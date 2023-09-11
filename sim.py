from pathlib import Path
from dataclasses import dataclass


@dataclass
class ColinearBlock:
    """
    dataclass to store a ColinearBlock for easier access.
    """

    genome1: str
    genome1_pos: tuple
    genome1_strand: str
    genome2: str
    genome2_pos: tuple
    genome2_strand: str


def read_xmfa(xmfa_file: Path) -> list[ColinearBlock]:
    """
    Reads xmfa entries only file and returns list of colinear blocks

    Args:
        xmfa_file (Path): Path to xmfa colinear block entries only file

    Returns:
        list[ColinearBlock]
    """
    with open(xmfa_file, "r") as f:
        lines = f.readlines()

    out = []
    for i, line in enumerate(lines):
        if line.strip().startswith("="):
            prev_two_lines = [lines[i - 2], lines[i - 1]]
            if all([prev_line.startswith(">") for prev_line in prev_two_lines]):
                prev_two_lines = [
                    prev_line.strip().split() for prev_line in prev_two_lines
                ]
                genome1_pos = tuple(
                    map(int, prev_two_lines[0][1].split(":")[1].split("-"))
                )
                genome2_pos = tuple(
                    map(int, prev_two_lines[1][1].split(":")[1].split("-"))
                )
                block = ColinearBlock(
                    genome1="wmel",
                    genome1_pos=genome1_pos,
                    genome1_strand=prev_two_lines[0][2],
                    genome2="wri",
                    genome2_pos=genome2_pos,
                    genome2_strand=prev_two_lines[1][2],
                )

                out.append(block)
    return out


if __name__ == "__main__":
    read_xmfa("./data/entries_only.xmfa")
