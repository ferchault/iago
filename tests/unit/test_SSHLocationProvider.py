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

    def test__connected(self):
        # create private key
        private_key = paramiko.RSAKey.generate(2048)

        # add to authorized keys
        authorized = os.path.join(os.path.expanduser('~'), '.ssh', 'authorized_keys')
        try:
            shutil.copyfile(authorized, authorized + '.bak')
        except IOError:
            pass # no authorized_keys file already present
        fh = open(authorized, 'a+')
        fh.write('\nssh-rsa ' + private_key.get_base64())
        fh.close()

        # open SSH connection to localhost with current user
        ssh = lp.SSHLocationProvider('localhost:' + os.path.join(BASE_DIR, 'fixtures'), pkey=private_key)

        # check presence of buckets
        bl = ssh.get_bucket_list()
        self.assertTrue('debug' in bl.name.values)

        # check file access
        self.assertTrue(ssh.has_file('debug-88ac57b3e437fa1d5d26d00d6c768324', 'index.ndx'))
        try:
            fh = ssh.open_file('debug-88ac57b3e437fa1d5d26d00d6c768324', 'index.ndx')
            fh.close()
        except:
            self.fail('Cannot open SSH remote file.')
        fn, emphemeral = ssh.local_file('debug-88ac57b3e437fa1d5d26d00d6c768324', 'index.ndx')
        self.assertTrue(emphemeral)
        os.remove(fn)

        # reset authorized keys
        try:
            shutil.copyfile(authorized + '.bak', authorized)
            os.remove(authorized + '.bak')
        except IOError:
            # no authorized_keys file present initially
            os.remove(authorized)
