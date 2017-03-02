from unittest import TestCase
from iago import LocationProvider as lp
import paramiko
import os

class TestSSHLocationProvider(TestCase):
    def test__prepare_connect(self):

        path = 'user@host:dir'
        client, port, sock, username, hostname, basepath = lp.SSHLocationProvider._prepare_connect(path)
        self.assertEqual(port, 22)
        self.assertEqual(username, 'user')
        self.assertEqual(hostname, 'host')
        self.assertEqual(basepath, 'dir')
        self.assertIsInstance(client, paramiko.SSHClient)

        path = 'foobar:dir'
        print os.getcwd()
        fh = open('.ssh_config', 'w')
        fh.write('Host foobar\n\tUser snafu\n\tHostName foobar.example.com\n')
        fh.close()
        client, port, sock, username, hostname, basepath = lp.SSHLocationProvider._prepare_connect(path)
        self.assertEqual(port, 22)
        self.assertEqual(username, 'snafu')
        self.assertEqual(hostname, 'foobar.example.com')
        self.assertEqual(basepath, 'dir')
        os.remove('.ssh_config')

        self.assertRaises(ValueError, lp.SSHLocationProvider._prepare_connect, 'foobar')

        self.assertRaises(ValueError, lp.SSHLocationProvider._prepare_connect, 'unknown_host:dir')


