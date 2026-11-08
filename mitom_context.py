"""mitom_context.py - RunContext: what the steps of one workflow run share (instead of module globals).

mitom_pipeline.main() makes one per command; resolve_tools (mitom_prepare) stores the --clade auto plan in it and the
workflows of mitom_run pass it to their steps. A function called without one gets a fresh RunContext, so a test, a
script or two runs in one process never see each other's state."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

from mitom_groups import Plan
from mitom_orient import OrientNotes


@dataclass
class RunContext:
    plan: Optional[Plan] = None                                   # --clade auto: phylum and code of every record and why
    dedup_stats: Dict[str, int] = field(default_factory=dict)     # GenBank input: records read / dropped (step 1 lines)
    orient: OrientNotes = field(default_factory=OrientNotes)      # fallback start genes, linear records reported
