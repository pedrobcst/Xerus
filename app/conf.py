from dataclasses import dataclass
from typing import Tuple


@dataclass
class AppSettings:
    TMP_FOLDER: str = "tmp"
    RESULTS_TMP_FOLDER: str = "results_tmp"
    KEEP_COLUMNS: Tuple = (
        "name",
        "rwp",
        "spacegroup",
        "crystal_system",
        "system_type"
    )

