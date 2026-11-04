"""How a run seed becomes a sub-simulation's seed.

A clonal run is not one simulation. Each clone has a founder, each
family may be re-grown several times before one survives, and each
observed cell or read is simulated on its own. Every one of those needs
a seed, and the run has one.

The first schedule derived them by addition: clone *c* took
``seed + c * 1_000_000``, a retry took ``clone_seed + attempt * 7919``,
a read or node took ``clone_seed + 1 + index``. Addition is a
translation, so that made the whole schedule a single lattice shared by
every run seed, with three consequences:

* **Runs reused each other's streams.** Clone *c* at run seed *S* took
  the seed clone 0 takes at run seed ``S + c * 1_000_000``, so the same
  founder -- and, with the read offsets agreeing too, the same corrupted
  read -- appeared in both runs. "Run it again at another seed" returned
  a reshuffling of clones already drawn rather than an independent
  replicate.
* **Two axes met.** A node took ``clone_seed + 1 + id`` and a retry took
  ``clone_seed + attempt * 7919``, so node 7918 and attempt 1 drew the
  same words.
* **A clone's block had no edge.** Node offsets grew with the tree and
  nothing stopped them crossing into the next clone's block.

The current schedule derives every axis through the engine's keyed
primitive instead, with the axis as the stream's role, so the run seed is
mixed rather than translated and the roles cannot meet.

Both schedules live here because the label is a wire value: a
family-replay bundle records which schedule drew it, and a bundle
written under the first one still replays.
"""

from typing import Dict, Optional


class _SeedSchedule:
    """The four seeds a clonal run derives, under one named scheme."""

    #: The name a family-replay bundle records, so a reader can pick the
    #: scheme its seeds were drawn under rather than assume the current one.
    label: str

    #: The scheme's own parameters that a bundle repeats as execution
    #: fields, beside the label. Only what was actually recorded: the
    #: additive scheme wrote down its retry stride and not its clone
    #: stride, and a replay has to accept the bundles that exist.
    recorded_parameters: Dict[str, int] = {}

    def unaddressable(self, run_seed: int, n_clones: int, retries: int) -> Optional[str]:
        """Why this scheme cannot seed that run, or ``None`` if it can."""
        return None

    def clone(self, run_seed: int, clone_index: int) -> int:
        """The seed for clone ``clone_index``'s founder."""
        raise NotImplementedError

    def attempt(self, clone_seed: int, attempt_index: int) -> int:
        """The seed for the ``attempt_index``-th try at growing a family."""
        raise NotImplementedError

    def read(self, clone_seed: int, read_index: int) -> int:
        """The seed for the ``read_index``-th draw off a founder.

        One axis for a repertoire's reads and for ``expand_clones``'
        descendants: both are a post-fork plan run against a founder,
        and no plan is both at once.
        """
        raise NotImplementedError

    def node(self, clone_seed: int, node_id: int) -> int:
        """The seed for observing the tree node ``node_id``."""
        raise NotImplementedError


class _Additive(_SeedSchedule):
    """The first schedule, kept so its bundles still replay.

    Reproduced exactly, including the shared ``+ 1 +`` offset that made
    a read and a node of the same index one stream: a replay has to
    re-derive the seeds the capture actually used, not better ones.
    """

    label = "clone_1000000_attempt_7919_node_plus_one_v1"
    recorded_parameters = {"retry_stride": 7_919}
    clone_stride = 1_000_000
    retry_stride = 7_919

    def unaddressable(self, run_seed, n_clones, retries):
        """Addition runs out of u64. The last clone's last retry is the
        highest seed the scheme reaches, and the engine takes a u64."""
        highest = self.attempt(self.clone(run_seed, n_clones - 1), retries)
        if highest > 2**64 - 1:
            return "clone/retry seed schedule exceeds u64"
        return None

    def clone(self, run_seed: int, clone_index: int) -> int:
        return int(run_seed) + clone_index * self.clone_stride

    def attempt(self, clone_seed: int, attempt_index: int) -> int:
        return clone_seed + attempt_index * self.retry_stride

    def read(self, clone_seed: int, read_index: int) -> int:
        return clone_seed + 1 + read_index

    def node(self, clone_seed: int, node_id: int) -> int:
        return clone_seed + 1 + node_id


class _Keyed(_SeedSchedule):
    """The current schedule: one keyed substream per axis.

    The axis name is the stream's role, so ``read`` 7918 and ``attempt``
    1 are different streams for the same reason that two passes of
    different names are, and a run seed is folded in rather than added
    to. ``node`` and ``read`` stay distinct roles even though no plan
    draws both: a tree's observation and a flat clone's read are
    different sub-simulations, and giving them one role would re-create
    the collision this schedule removes the moment a plan does both.
    """

    label = "keyed_substream_clone_attempt_read_node_v2"

    def unaddressable(self, run_seed, n_clones, retries):
        """Nothing runs out: a derived seed is a whole u64 whatever went
        in. The one limit is how far an axis counts, since a stream's
        role carries a 32-bit occurrence."""
        if n_clones > 2**32:
            return "more clones than one stream role can name"
        return None

    @staticmethod
    def _derive(seed: int, role: str, index: int) -> int:
        import GenAIRR._engine as _engine

        return _engine.substream_seed(int(seed), role, index)

    def clone(self, run_seed: int, clone_index: int) -> int:
        return self._derive(run_seed, "clone", clone_index)

    def attempt(self, clone_seed: int, attempt_index: int) -> int:
        return self._derive(clone_seed, "attempt", attempt_index)

    def read(self, clone_seed: int, read_index: int) -> int:
        return self._derive(clone_seed, "read", read_index)

    def node(self, clone_seed: int, node_id: int) -> int:
        return self._derive(clone_seed, "node", node_id)


#: The schedule every new run draws under.
CURRENT = _Keyed()

#: Every schedule a recorded label can name, for replaying an artifact
#: that was not drawn under :data:`CURRENT`.
BY_LABEL: Dict[str, _SeedSchedule] = {
    _Additive.label: _Additive(),
    _Keyed.label: CURRENT,
}
