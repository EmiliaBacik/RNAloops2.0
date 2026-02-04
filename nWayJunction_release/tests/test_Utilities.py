from unittest import TestCase

from Commons.Utilities import pairs_generator

class TestUtilities(TestCase):

    def test_pairs_generator(self):
        list = [1,2,3,4]
        test_list = pairs_generator(list)
        self.assertEqual(tuple(test_list), ((1, 2), (2, 3), (3, 4), (4, 1)))
