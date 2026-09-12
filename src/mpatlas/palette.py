"""Atlas figure color palette."""

from __future__ import annotations

PALETTE = {
    "sage": "#8FB996",
    "teal": "#5FA8A8",
    "peach": "#F4A261",
    "coral": "#E76F51",
    "lavender": "#B8A9E8",
    "mustard": "#E9C46A",
    "navy": "#264653",
}

TEXT = PALETTE["navy"]
MUTED = "#5A6F74"
GRID = "#D9E2E0"
FACE = "#FFFFFF"
SPINE = "#A8B8B5"
POS = PALETTE["sage"]
NEG = PALETTE["coral"]
MID = PALETTE["mustard"]

SPLIT_COLORS = {
    "random": PALETTE["sage"],
    "cluster": PALETTE["coral"],
    "center": PALETTE["peach"],
    "organism": PALETTE["teal"],
    "architecture": PALETTE["lavender"],
    "time": PALETTE["mustard"],
}

SERIES = [
    PALETTE["teal"],
    PALETTE["sage"],
    PALETTE["mustard"],
    PALETTE["peach"],
    PALETTE["coral"],
    PALETTE["lavender"],
]
