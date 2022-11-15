import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd


def read_metadata(file_name: Path):
    return pd.read_csv(file_name, sep="\t", header=None)


def clear_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    def parse_timepoint(string: str) -> Optional[str]:
        if "T0" in string:
            return "0"
        elif "T6" in string:
            return "6"
        else:
            return None

    def parse_group(string: str) -> Optional[str]:
        if "NS" in string:
            return "Non-survivor"
        elif re.match("S[0-9]{1}", string):
            return "Survivor"
        elif "HC" in string:
            return "Healthy"

    metadata.columns = pd.Index(["GEO", "Group"])
    metadata.loc[:, "timepoint"] = metadata["Group"].apply(parse_timepoint)
    metadata.loc[:, "group"] = metadata["Group"].apply(parse_group)
    return metadata


def main(): 
    metadata_path = Path(sys.argv[1])
    if metadata_path is None:
        metadata_path = Path("data/raw/metadata.txt")
    metadata = read_metadata(metadata_path)
    metadata = clear_metadata(metadata)
    metadata.to_csv("data/processed/metadata.csv", index=False)


if __name__ == "__main__":
    main()
