"""Compare the local Cl V source with an unverified indexed journal table.

The table transcription is from a search index of the typeset 1999 article.
It is diagnostic evidence only until the final article is read directly.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


SOURCE = Path(
    "/home/kiseon/RT_Codes/MOCASSIN/"
    "mocassin-mocassin.2.02.73.2/data/clv.dat"
)
INDEXED_TABLE = {
    1000: 2.59,
    1500: 2.63,
    2000: 2.72,
    2500: 2.83,
    3000: 2.94,
    3500: 3.05,
    4000: 3.15,
    5000: 3.33,
    6000: 3.49,
    7000: 3.62,
    8000: 3.73,
    9000: 3.82,
    10000: 3.90,
    11000: 3.98,
    12000: 4.04,
    13000: 4.10,
    14000: 4.16,
    15000: 4.21,
    16000: 4.26,
    18000: 4.35,
    20000: 4.43,
    24000: 4.57,
    28000: 4.68,
    32000: 4.78,
    36000: 4.85,
}
INDEX_URL = (
    "https://citeseerx.ist.psu.edu/document?"
    "doi=15e037d5725168fa42aab91856e8d82d4ab83d20"
    "&repid=rep1&type=pdf"
)


def source_rows() -> list[dict]:
    lines = SOURCE.read_text().splitlines()
    assert int(lines[0]) == 5
    assert [int(v) for v in lines[6].split()] == [2, 15, 0]
    temperatures = [int(v) for v in lines[9:24]]
    strengths = [float(line.split()[2]) for line in lines[25:40]]
    assert len(temperatures) == len(strengths) == 15
    assert temperatures == sorted(set(temperatures))
    return [dict(temperature_K=t, upsilon=u) for t, u in zip(temperatures, strengths)]


def main() -> None:
    rows = source_rows()
    comparison = [
        {
            "temperature_K": row["temperature_K"],
            "mocassin_upsilon": row["upsilon"],
            "indexed_table_upsilon": INDEXED_TABLE[row["temperature_K"]],
            "absolute_difference": abs(
                row["upsilon"] - INDEXED_TABLE[row["temperature_K"]]
            ),
        }
        for row in rows
        if row["temperature_K"] in INDEXED_TABLE
    ]
    result = {
        "source_path": str(SOURCE),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_nodes": rows,
        "source_node_count": len(rows),
        "source_temperature_range_K": [rows[0]["temperature_K"], rows[-1]["temperature_K"]],
        "indexed_journal_table_url": INDEX_URL,
        "indexed_table_direct_final_pdf_verified": False,
        "indexed_table_transcribed_row_count": len(INDEXED_TABLE),
        "indexed_table_last_40000K_row_incomplete": True,
        "exact_temperature_comparison": comparison,
        "exact_temperature_overlap_count": len(comparison),
        "maximum_absolute_overlap_difference": max(
            row["absolute_difference"] for row in comparison
        ),
        "cold_constant_1000K": rows[0]["upsilon"],
        "indexed_table_1000K": INDEXED_TABLE[1000],
        "cold_constant_1000K_fractional_difference": (
            rows[0]["upsilon"] / INDEXED_TABLE[1000] - 1
        ),
        "publication_source_accepted": False,
    }
    output = Path(__file__).with_name("result.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(output)
    print("overlap", len(comparison), "max_abs", result["maximum_absolute_overlap_difference"])
    print("1000K cold fractional difference", result["cold_constant_1000K_fractional_difference"])


if __name__ == "__main__":
    main()
