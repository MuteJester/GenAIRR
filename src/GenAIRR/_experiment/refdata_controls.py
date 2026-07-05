from __future__ import annotations


class _RefdataControlsMixin:
    __slots__ = ()

    def curate_refdata(
        self,
        policy: str,
        *,
        allowed=None,
        keep_unannotated: bool = True,
    ) -> "Experiment":
        """**Curation** — select which subset of the catalogue
        participates in simulation. Replaces this experiment's
        reference cartridge with a curated version.

        In the cartridge model, curation is distinct from validation
        and from the catalogue itself:

        - **validation** describes what the catalogue contains;
        - **curation** decides which alleles participate;
        - **simulation** runs against the curated cartridge.

        ``policy`` is one of:

        - ``"raw"`` — identity policy; no-op (kept for symmetry).
        - ``"functional_anchors_only"`` — drop V and J alleles whose
          anchor doesn't satisfy the active anchor rule (missing
          anchor, codon out of bounds, or codon AA outside the
          rule's ``expected_amino_acids``). D and C pools pass
          through unchanged.
        - ``"functional_status"`` — filter V/D/J pools by IMGT
          functional status. ``allowed`` is the list of statuses
          to keep (default ``["functional"]``); accepted strings
          are ``"functional"``, ``"orf"``, ``"pseudogene"``, and
          ``"unknown"`` (case-insensitive, ``"F"`` / ``"P"``
          aliases also work). ``keep_unannotated`` (default
          ``True``) controls whether alleles with no annotation
          survive — bundled ``.pkl`` cartridges currently leave
          status unannotated, so the default preserves backward
          compatibility. C pool passes through unchanged.

        Use this when the catalogue (e.g. bundled ``mouse_igh`` or
        ``human_tcrb``) includes pseudogene/ORF alleles you'd rather
        exclude than permit at runtime via
        :meth:`allow_curatable_refdata`. Curation is the
        professional model; ``allow_curatable_refdata`` is the
        broader runtime opt-in for sampling the raw catalogue
        as-is.

        Curation cannot fix structural problems — duplicate allele
        names, invalid sequence bytes, and locus/chain-type
        mismatches still surface from the compile-time validator
        regardless of the curated cartridge state. If curation
        empties a required pool, ``compile()`` fails with
        ``EmptyRequiredPool``.

        The curated cartridge's ``identity.source`` is tagged
        ``|curated:<policy>`` so trace files and content hashes
        distinguish raw from curated artefacts. Returns ``self`` so
        the call chains fluently. See ``docs/reference_cartridge.md``.
        """
        kwargs = {}
        if allowed is not None:
            kwargs["allowed"] = list(allowed)
        if policy == "functional_status":
            kwargs["keep_unannotated"] = bool(keep_unannotated)
        self._refdata = self._refdata.curated(policy, **kwargs)
        return self

    def allow_curatable_refdata(self, enabled: bool = True) -> "Experiment":
        """Opt in to the lenient `AllowCuratable` refdata validation
        mode for subsequent ``compile`` / ``run`` / ``run_records``
        calls. Returns ``self`` so the call chains fluently.

        Sits alongside :meth:`curate_refdata` as the cartridge's two
        ways to handle pseudogene-bearing catalogues:

        - ``curate_refdata("functional_anchors_only")`` **removes**
          the non-canonical alleles. Strict validation then passes
          because the curated catalogue is clean. This is the
          professional model — the cartridge identifies which
          alleles actually participate.
        - ``allow_curatable_refdata()`` **keeps** the catalogue
          as-is and relaxes the validator. Strict validation passes
          Curatable issues (pseudogene-shape anchor anomalies) but
          still rejects Fatal ones (empty pools, duplicates, invalid
          bytes, anchor out of bounds, locus/chain mismatch).

        Fatal issues are never opt-outable. Curatable issues — V
        anchor codon not Cys, J anchor codon outside the locus's
        expected set, missing V/J anchor — reflect pseudogene/ORF
        entries in real reference catalogues (the bundled
        ``mouse_igh`` and ``human_tcrb`` data both contain them).

        Recommended progression: start strict; if your catalogue
        contains pseudogenes you want to *exclude*, use
        :meth:`curate_refdata`; if you want to *sample from* them
        explicitly, use this method. See
        ``docs/reference_cartridge.md``.
        """
        self._allow_curatable_refdata = bool(enabled)
        return self
