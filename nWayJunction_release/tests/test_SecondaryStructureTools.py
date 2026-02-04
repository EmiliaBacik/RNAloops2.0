from unittest import TestCase

from Commons.SecondaryStructureTools import find_matching_parenthesis, find_junction, find_matching_char, fill_secondary

class TestSecondaryStructureTools(TestCase):
    def test_find_matching_parenthesis_1(self):
        test_string = '(........)'
        self.assertEqual(find_matching_parenthesis(test_string, 0), 9)

    def test_find_matching_parenthesis_2(self):
        test_string = '(....(.)...)'
        self.assertEqual(find_matching_parenthesis(test_string, 0), 11)

    def test_find_matching_parenthesis_3(self):
        test_string = '(....(.)...)'
        self.assertEqual(find_matching_parenthesis(test_string, 5), 7)

    def test_find_junction_1(self):
        test_string = '(..((.)).(.)...)'
        self.assertTrue(find_junction(test_string, 0)[0])
        self.assertEqual(find_junction(test_string, 0)[1], 3)
        self.assertEqual(find_junction(test_string, 0)[2], [[1, 16], [4, 8], [10, 12]])

    def test_find_junction_2(self):
        test_string = '..........'
        self.assertFalse(find_junction(test_string, 0)[0])
        self.assertEqual(find_junction(test_string, 0)[1], 0)
        self.assertEqual(find_junction(test_string, 0)[1], 0)

    def test_find_matching_char_1(self):
        test_string = '(.[.((.])).(.)...)'
        self.assertEqual(find_matching_char('[', test_string, 2 ),7)
        self.assertEqual(find_matching_char(']', test_string, 7 ),2)

    def test_fill_secondary_1(self):
        test_string = '(.[.((.])).(.)...)'
        connectors = [[0,3],[4,4],[8,10]]
        self.assertEqual(fill_secondary(connectors, test_string)[0], ['(..','','))'])
        self.assertEqual(fill_secondary(connectors, test_string)[1], '(...((..)).(.)...)')

    def test_fill_secondary_2(self):
        test_string = '(.[.((.])).(.)...)'
        connectors = [[0,3],[4,8],[9,10]]
        self.assertEqual(fill_secondary(connectors, test_string)[0], ['(.[', '((.]', ')'])
        self.assertEqual(fill_secondary(connectors, test_string)[1], '(.[.((.])).(.)...)')


