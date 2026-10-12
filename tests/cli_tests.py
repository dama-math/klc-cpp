"""Check CLI error reporting and successful output without third-party packages."""
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    binary = Path(sys.argv[1]).resolve()

    def run(*args):
        return subprocess.run([str(binary), *map(str, args)], capture_output=True,
                              text=True, timeout=30)

    for args in [("--random", "-1"), ("--random", "25"),
                 ("--random", "8", "-1"), ("--speedtest2", "12abc"),
                 ("--speedtest2", "999999999999999999999"),
                 ("--random", "8", "1", "42", "extra"),
                 ("--cyclic", "-1"), ("--unknown", "x"),
                 ("--file", "missing-file-4edb7.txt"), ("0123", "01")]:
        result = run(*args)
        assert result.returncode != 0 and result.stderr, (args, result)
    assert run("0123", "2301").stdout == "0123 2301 1,1\n"
    assert run("--simple", "2301").stdout == "2\n"
    with tempfile.TemporaryDirectory(dir=binary.parent) as tmp:
        path = Path(tmp) / "verify.txt"
        path.write_text("0123 2301 1,2\n", encoding="utf-8")
        result = run("--verify", path)
        assert result.returncode != 0 and "[FAIL]" in result.stdout
        path.write_text("# comment\n\n0123 2301 1,1\n", encoding="utf-8")
        result = run("--verify", path)
        assert result.returncode == 0 and "1/1 passed" in result.stdout
        path.write_text("malformed\n0123 2301 1,1\n", encoding="utf-8")
        result = run("--verify", path)
        assert result.returncode != 0 and "Error processing line" in result.stderr
        path.write_text("0123 2301\ninvalid 01\n", encoding="utf-8")
        result = run("--file", path)
        assert result.returncode != 0 and result.stdout == "0123 2301 1,1\n"
    print("CLI tests passed.")


if __name__ == "__main__":
    main()
