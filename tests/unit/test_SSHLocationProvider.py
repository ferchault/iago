from unittest import TestCase
from iago import LocationProvider as lp
import paramiko
import shutil
import os

BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


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
        fh = open('.ssh_config', 'w')
        fh.write('Host foobar\n\tUser snafu\n\tPort 24\n\tHostName foobar.example.com\n')
        fh.close()
        client, port, sock, username, hostname, basepath = lp.SSHLocationProvider._prepare_connect(path)
        self.assertEqual(port, 24)
        self.assertEqual(username, 'snafu')
        self.assertEqual(hostname, 'foobar.example.com')
        self.assertEqual(basepath, 'dir')
        os.remove('.ssh_config')

        self.assertRaises(ValueError, lp.SSHLocationProvider._prepare_connect, 'foobar')

    def test__connect(self):
        # create private key
        private_key = paramiko.RSAKey.generate(2048)
        key_filename = os.path.join(os.path.expanduser('~'), 'iago-sshkey')
        private_key.write_private_key_file(key_filename)
        fh = open(key_filename + '.pub', 'w').write(private_key.get_base64())

        # add to authorized keys
        authorized = os.path.join(os.path.expanduser('~'), '.ssh', 'authorized_keys')
        shutil.copyfile(authorized, authorized + '.bak')
        fh = open(authorized, 'a')
        fh.write('\n' + private_key.get_base64())
        fh.close()

        # open SSH connection to localhost with current user
        ssh = lp.SSHLocationProvider('localhost:' + os.path.join(BASE_DIR, 'fixtures'))
        bl = ssh.get_bucket_list()
        self.assertTrue('debug' in bl.name.values)

        # reset authorized keys
        shutil.copyfile(authorized + '.bak', authorized)
        os.remove(authorized + '.bak')
