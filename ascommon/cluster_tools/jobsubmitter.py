""".

TODO: Use a SQLite database if you ever want to run on SciNet.
"""
import os
import os.path as op
import time
import logging
import paramiko

logger = logging.getLogger(__name__)


class JobSubmitter:
    """.

    .. note::

        While it would be safer call the Python executable directly,
        rather than calling bash, this does not seem to be supported on PBS
        (at least not on banting :()).
    """

    def __init__(
        self, job_name, head_node_type, head_node_ip, log_root_path=None, *,
        concurrent_job_limit=None,
        queue='medium',
        nproc=1,
        walltime='24:00:00',
        mem='2000M',
        email='noname@example.com',
        email_opts='ae',
        qsub_shell='/bin/bash',  # Python does not work on Banting
        qsub_script=op.join(op.dirname(op.abspath(__file__)), 'scripts', 'qsub.sh'),
        qsub_script_args=None
    ):
        """.

        Parameters
        ----------
        job_name : str
            Name of the job.
        head_node_ip : str
            IP address of the head node.
        head_node_type : str
            Job manager type ['sge', 'pbs', ...]
        log_root_path : str, default None
            Location where the log files should be saved.
            A new folder will be created here for each job.
            TODO: Allow this to be a database.
        """
        # Required arguments
        self.job_name = job_name
        self.head_node_ip = head_node_ip

        if head_node_type not in ['sge', 'pbs']:
            raise ValueError("Wrong 'head_node_type': '{}'".format(head_node_type))
        self.head_node_type = head_node_type

        if log_root_path is None:
            log_root_path = op.join(op.expanduser('~'), 'pbs-output')
        else:
            log_root_path = op.abspath(log_root_path)
        os.makedirs(log_root_path, exist_ok=True)

        log_path = op.join(log_root_path, job_name)
        try:
            os.mkdir(log_path)
        except FileExistsError:
            logger.error("Each job needs to create it's own empty log folder!")
            raise
        time.sleep(1)
        self.log_path = log_path

        self.ssh = None

        # Default arguments
        self.concurrent_job_limit = concurrent_job_limit
        self.queue = queue
        self.nproc = nproc
        self.walltime = walltime
        self.mem = mem
        self.email = email
        self.email_opts = email_opts
        self.qsub_shell = qsub_shell  # Python does not work on Banting
        self.qsub_script = qsub_script
        self.qsub_script_args = qsub_script_args  # key-value dict

    # === Manage connection ===

    def connect(self):
        """Open connection to head node."""
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(self.head_node_ip)
        self.ssh = ssh

    def close(self):
        """Close connection to head node."""
        self.ssh.close()

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    # === QSUB ===

    @property
    def qsub_system_command_template(self):
        """Template of the system command which will submit the given job on the head node."""
        # ## Other possible options
        # -l s_vmem=5650M
        # -l h_vmem=5850M
        # -l mem_free=5850M
        # -l virtual_free=5850M
        # -cwd
        # -p 0  # priority
        # -v PATH=/home/kimlab1/strokach/anaconda3/bin
        _template = """\
qsub -S "{qsub_shell}" -N "{job_name}" -M "{email}" -m {email_opts} -V \
{_qsub_options_str} "{qsub_script}" {_qsub_script_args_str} \
"""
        template = _template.format(**{
            '_qsub_options_str': self._qsub_options_str,
            '_qsub_script_args_str': self._qsub_script_args_str,
            **self.__dict__
        })
        return template

    @property
    def _qsub_options_str(self):
        if self.head_node_type == 'sge':
            return """\
-pe smp {nproc} -q {queue} -l s_rt={walltime},h_rt={walltime},mem_free={mem} \
-v SYSTEM_COMMAND="{{system_command}}" \
-v STDOUT_LOG="{log_path}/{{job_id}}.out" \
-v STDERR_LOG="{log_path}/{{job_id}}.err" \
""".format(**self.__dict__)
        elif self.head_node_type == 'pbs':
            return """\
-l nodes=1:ppn={nproc},walltime={walltime},mem={mem} \
-v SYSTEM_COMMAND="{{system_command}}",\
STDOUT_LOG="{log_path}/{{job_id}}.out",\
STDERR_LOG="{log_path}/{{job_id}}.err" \
""".format(**self.__dict__)

    @property
    def _qsub_script_args_str(self):
        if self.qsub_script_args is None:
            return ''
        script_args_str = ' '.join(
            '--{} {}'.format(key, value) for key, value in self.qsub_script_args.items()
        )
        if self.head_node_type == 'sge':
            return script_args_str
        elif self.head_node_type == 'pbs':
            return '-F "{}"'.format(script_args_str)

    # === Submit jobs ===

    def submit(self, iterable, **kwargs):
        results = []
        for i, (job_id, system_command) in enumerate(iterable):
            self._respect_concurrent_job_limit(i)
            system_command = self.qsub_system_command_template.format(
                job_id=job_id,
                system_command=system_command.replace('"', '\\"'),
            )
            logger.debug(system_command)
            stdout, stderr = self._exec_system_command(system_command)
            time.sleep(0.05)
            results.append((stdout, stderr, ))
        return results

    def get_num_running_jobs(self):
        system_command = 'qstat -u "$USER" | grep "$USER" | grep -i " r  " | wc -l'
        stdout, stderr = self._exec_system_command(system_command)
        num_running_jobs = int(stdout)
        return num_running_jobs

    def _exec_system_command(self, system_command):
        n_tries = 0
        stdout = 'x.x'
        stderr = 'x.x'
        while n_tries < 5 and stderr:
            if n_tries:
                delay = (1 * n_tries)
                logger.debug('Sleeping for {} seconds...'.format(delay))
                time.sleep(delay)
            n_tries += 1
            stdin_fh, stdout_fh, stderr_fh = self.ssh.exec_command(system_command)
            stdout = stdout_fh.read().decode().strip()
            if stdout:
                logger.debug(stdout)
            stderr = stderr_fh.read().decode().strip()
            if stderr:
                logger.error(stderr)
        return stdout, stderr

    def _respect_concurrent_job_limit(self, i):
            # Limit the number of jobs running simultaneously
            STEP_SIZE = 50
            DELAY = 120
            if self.concurrent_job_limit is not None and ((i + 1) % STEP_SIZE) == 0:
                num_running_jobs = self.get_num_running_jobs()
                while (num_running_jobs + STEP_SIZE) > self.concurrent_job_limit:
                    logger.info(
                        "'num_running_jobs' reached! Sleeping for {:.0f} minutes..."
                        .format(DELAY / 60))
                    time.sleep(DELAY)
                    num_running_jobs = self.get_num_running_jobs()
