from pathlib import Path
from dataclasses import dataclass


@dataclass
class ColinearBlock:
    """
    dataclass to store a ColinearBlock for easier access.
    """

    genome1: str
    genome1_pos: str
    genome1_strand: str
    genome2: str
    genome2_pos: str
    genome2_strand: str


def read_xmfa(xmfa_file: Path) -> list[ColinearBlock]:
    """
    Reads xmfa entries only file and returns list of tuples of colinear blocks

    Args:
        xmfa_file (Path): Path to xmfa colinear block entries only file

    Returns:
        list[tuple[tuple]]: List of tuples of colinear blocks
                        ex) [((wRi,1,100), (wmel,5,100))...]
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
                block = ColinearBlock(
                    genome1="wmel",
                    genome1_pos=prev_two_lines[0][1],
                    genome1_strand=prev_two_lines[0][2],
                    genome2="wri",
                    genome2_pos=prev_two_lines[1][1],
                    genome2_strand=prev_two_lines[1][2],
                )
                out.append(block)
    return out


if __name__ == "__main__":
    read_xmfa("./data/entries_only.xmfa")
