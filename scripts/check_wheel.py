"""Build a wheel and test the *installed* package, not the source tree.

Running the suite from a checkout proves the source works. It does not prove
the wheel does: a module left out of the package, a missing entry point or a
data file that never made it into the distribution all pass a source-tree run
and fail for whoever installs it.

This builds the distributions, installs the wheel into a throwaway virtual
environment with ``PYTHONPATH`` stripped so nothing resolves back to ``src/``,
and runs the full suite plus a command-line smoke test from a directory that is
not the repository.

    python scripts/check_wheel.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run(
    command: list[str], cwd: Path, env: dict[str, str]
) -> subprocess.CompletedProcess:
    print(f"$ {' '.join(str(part) for part in command)}")
    return subprocess.run(command, cwd=cwd, env=env, text=True)


def main() -> int:
    clean = {key: value for key, value in os.environ.items() if key != "PYTHONPATH"}

    dist = ROOT / "dist"
    if dist.exists():
        shutil.rmtree(dist)
    built = run([sys.executable, "-m", "build", "--no-isolation"], ROOT, clean)
    if built.returncode:
        print("build failed", file=sys.stderr)
        return 1

    wheels = sorted(dist.glob("*.whl"))
    if len(wheels) != 1:
        print(f"expected exactly one wheel, found {len(wheels)}", file=sys.stderr)
        return 1
    wheel = wheels[0]

    with tempfile.TemporaryDirectory(prefix="plasmid-map-wheel-") as holding:
        home = Path(holding)
        environment = home / "venv"
        venv.create(environment, with_pip=True)
        binaries = environment / ("Scripts" if os.name == "nt" else "bin")
        python = binaries / "python"

        install = run(
            [str(python), "-m", "pip", "install", "--quiet", f"{wheel}[dev]"],
            home,
            clean,
        )
        if install.returncode:
            print("install failed", file=sys.stderr)
            return 1

        # Tests and examples are not installed by the wheel, so they are run
        # from the checkout while the package itself resolves to the venv.
        suite = run(
            [str(python), "-m", "pytest", "-q", str(ROOT / "tests")], ROOT, clean
        )
        if suite.returncode:
            print("installed package fails its own test suite", file=sys.stderr)
            return 1

        examples = ROOT / "examples"
        command = str(binaries / "plasmid-map")
        plasmid = str(examples / "pbr322.gb")
        phage = str(examples / "phix174.gb")
        virus = str(examples / "aav2.gb")
        cassette = str(examples / "cassette.tsv")
        checks: list[list[str]] = [
            [command, "--version"],
            [command, "types"],
            [command, "features", plasmid],
            [command, "circular", plasmid],
            [command, "circular", plasmid, "-f", "json"],
            [
                command,
                "circular",
                plasmid,
                "--kind",
                "CDS,rep_origin,repeat_region,regulatory",
                "--gc",
                "--legend",
                "--width",
                "900",
                "--svg",
                str(home / "pbr322.svg"),
                "--output",
                str(home / "pbr322.tsv"),
                "--report",
                str(home / "pbr322.html"),
                "-f",
                "none",
            ],
            [
                command,
                "circular",
                phage,
                "--kind",
                "CDS",
                "--start",
                "3981",
                "--svg",
                str(home / "phix174.svg"),
                "-f",
                "none",
            ],
            [
                command,
                "linear",
                virus,
                "--kind",
                "CDS,repeat_region",
                "--gc",
                "--width",
                "900",
                "--svg",
                str(home / "aav2.svg"),
                "-f",
                "none",
            ],
            [
                command,
                "linear",
                "--length",
                "2744",
                "--features",
                cassette,
                "--title",
                "AAV-CMV-EGFP cassette",
                "--svg",
                str(home / "cassette.svg"),
                "-f",
                "none",
            ],
            [str(python), "-m", "plasmid_map", "--version"],
        ]
        for check in checks:
            smoke = run(check, home, clean)
            if smoke.returncode:
                print(f"smoke test failed: {' '.join(check)}", file=sys.stderr)
                return 1

        for name in (
            "pbr322.svg",
            "pbr322.tsv",
            "pbr322.html",
            "phix174.svg",
            "aav2.svg",
            "cassette.svg",
        ):
            if not (home / name).exists():
                print(f"the run did not write {name}", file=sys.stderr)
                return 1
            if name.endswith(".svg"):
                body = (home / name).read_text(encoding="utf-8")
                if not body.startswith("<?xml") or "<svg" not in body:
                    print(f"{name} is not an SVG document", file=sys.stderr)
                    return 1

    print("\nwheel verified: installed package passes the suite and the CLI runs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
