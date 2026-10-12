"""Check fixture generation failure handling and reproducibility."""
import importlib.util
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("generate_verify", ROOT / "generate_verify.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)


class GeneratorTests(unittest.TestCase):
    def test_deterministic_quotas(self):
        def calculator(binary, x, w, timeout):
            return ("0", "1", "1,1")[int(x[0]) % 3]
        group = generator.GROUPS[0]
        with patch.object(generator, "run_klc", side_effect=calculator):
            first = generator.generate_group(None, random.Random(42), group, 1000, 1)
            second = generator.generate_group(None, random.Random(42), group, 1000, 1)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 16)
        coefficients = [line.split()[2] for line in first[1:]]
        self.assertEqual([coefficients.count(p) for p in ("0", "1", "1,1")], [2, 5, 8])

    def test_failure_preserves_existing_output(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            output = Path(tmp) / "fixture.txt"
            output.write_text("previous fixture", encoding="utf-8")
            with patch.object(generator, "run_klc", return_value=None):
                with self.assertRaises(RuntimeError):
                    generator.generate(None, output, 42, 1, 2, [generator.GROUPS[0]])
            self.assertEqual(output.read_text(encoding="utf-8"), "previous fixture")

    def test_notation_and_pattern(self):
        self.assertEqual(generator.to_one_line(range(12)), "0123456789ab")
        self.assertFalse(generator.is_3412_avoiding([2, 3, 0, 1]))
        self.assertTrue(generator.is_3412_avoiding([0, 1, 2, 3]))


if __name__ == "__main__":
    unittest.main()
