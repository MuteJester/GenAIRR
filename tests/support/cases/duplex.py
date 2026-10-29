"""Short authored controls shared by duplex Experiment integration tests."""
from GenAIRR import RefDataConfig
from GenAIRR.experimental import CoupledDuplexRepair, DuplexJoin, DuplexHistory, DuplexRepairModel


def references(v_count=1):
    r = RefDataConfig.vdj()
    for i in range(v_count):
        r.add_v_allele(f"V{i+1}*01", f"V{i+1}", b"AAATGT", anchor=3)
    r.add_d_allele("D*01", "D", b"CCC")
    r.add_j_allele("J*01", "J", b"TTTAAA", anchor=0)
    return r


def script(*, inverted=False, order=("vd", "dj"), restore=False, added=""):
    s = CoupledDuplexRepair("AAATGT", "CCC", "TTTAAA", d_inverted=inverted)
    for name in order:
        join = DuplexJoin(name).open("left", offset=0).open("right", offset=0)
        if name == "vd":
            if restore:
                join = join.resect("left", "top", count=1)
            if added:
                mate = added.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                join = join.add_n("left", added).add_n("right", mate)
            join = join.synapse(right_start=6 + len(added))
            if restore:
                join = join.fill("top")
        else:
            join = join.synapse(right_start=3)
        s = s.join(join.ligate())
    return s


def history(name="blunt", *, weight=1.0, v_call="V1*01", **kwargs):
    return DuplexHistory(name, script(**kwargs), v_call, "D*01", "J*01", weight)


def model(ref, rows=None, **kwargs):
    return DuplexRepairModel.bind(ref, rows or [history()], reference_ends_are_coding_ends=True, **kwargs)
