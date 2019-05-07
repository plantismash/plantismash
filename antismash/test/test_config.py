try:
    import unittest2
except ImportError:
    import unittest as unittest2
from os import path
from argparse import Namespace
from antismash import config


class TestConfig(unittest2.TestCase):
    def setUp(self):
        self.original_config = config._config
        self.original_basedir = config._basedir
        self.original_default_name = config._default_name
        self.original_sys_name = config._sys_name
        self.original_user_file_name = config._user_file_name

    def tearDown(self):
        config._config = self.original_config
        config._basedir = self.original_basedir
        config._default_name = self.original_default_name
        config._sys_name = self.original_sys_name
        config._user_file_name = self.original_user_file_name

    def test_load_config(self):
        "Test config.load_config()"
        config._basedir = path.dirname(__file__)
        config._default_name = 'test.cfg'
        config._sys_name = 'test_sys.cfg'
        config._user_file_name = path.join(config._basedir, 'test_user.cfg')
        c = Namespace(testing=True)
        config.load_config(c)
        self.assertTrue(c.testing)
        self.assertEqual('true', c.default_file_loaded)
        self.assertEqual('true', c.sublevel.exists)
        self.assertEqual('true', c.sys_file_loaded)
        self.assertEqual('true', c.user_file_loaded)

    def test_set_config(self):
        "Test config.set_config()"
        c = Namespace(testing=True)
        self.assertIsNone(config._config)
        config.set_config(c)
        self.assertEqual(c, config._config)

    def test_get_config(self):
        "Test config.get_config()"
        config._config = Namespace(testing=True)
        c = config.get_config()
        self.assertEqual(config._config, c)
