import unittest

from fsac import update

class TestUpdate(unittest.TestCase):

    def test_reverse_complement(self):

        correct_answer = 'NTGTAATC'
        result = update.reverse_complement('GATTACAN')

        self.assertEqual(result, correct_answer)

def main():
    unittest.main()


if __name__ == '__main__':
    main()
