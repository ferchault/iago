from unittest import TestCase
from iago import LocationProvider as lp
import paramiko

class TestSSHLocationProvider(TestCase):
    def test__prepare_connect(self):

        path = 'user@host:dir'
        client, port, sock, username, hostname, basepath = lp.SSHLocationProvider._prepare_connect(path)
        self.assertEqual(port, 22)
        self.assertEqual(username, 'user')
        self.assertEqual(hostname, 'host')
        self.assertEqual(basepath, 'dir')

        self.assertIsInstance(client, paramiko.SSHClient)