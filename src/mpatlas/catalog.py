from __future__ import annotations

from dataclasses import asdict, dataclass, field

QUESTIONS = ("A", "B", "B-conditions", "C")

STATUS_RANK = {
    "selected": 1,
    "cloned": 2,
    "expressed": 3,
    "soluble": 4,
    "purified": 5,
    "crystallized": 6,
    "diffraction": 6,
    "nmr": 6,
    "in pdb": 7,
    "structures": 7,
    "in_pdb": 7,
}


def normalize_status(raw: str | None) -> str | None:
    if not raw:
        return None
    s = raw.strip().lower().replace("_", " ").replace("-", " ")
    aliases = {
        "work stopped": "selected",
        "target selected": "selected",
        "pdb": "in pdb",
        "in pdb": "in pdb",
        "structure determined": "in pdb",
        "crystal structure": "crystallized",
        "diffraction-quality crystals": "diffraction",
        "diffraction quality crystals": "diffraction",
        "hsqc": "nmr",
        "nmr assigned": "nmr",
        "nmr structure": "nmr",
        "membrane protein solubilized": "soluble",
        "expression tested": "expressed",
    }
    s = aliases.get(s, s)
    return s if s in STATUS_RANK else s


@dataclass
class Record:
    source: str
    question: str
    sequence: str
    label: float | None = None
    status: str | None = None
    organism: str | None = None
    center: str | None = None
    accession: str | None = None
    pdb_id: str | None = None
    architecture: str | None = None
    n_tm: int | None = None
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        extra = d.pop("extra") or {}
        d.update({k: v for k, v in extra.items() if k not in d})
        return d
