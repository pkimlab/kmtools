# -*- coding: utf-8 -*-
import subprocess
import six
import time

# %%

def submit_job(system_command):
    """Submit a job by running ``system_command``.
    """
    n_tries = 0
    jobnumber = 0
#    print(system_command)
    if cluster_type is None and system_command.strip().startswith('submitjob'):
        system_command = 'echo ' + system_command
    while not jobnumber and n_tries < 720:
        n_tries += 1
        time.sleep(0.5) # BlockingIOError: [Errno 11] Resource temporarily unavailable
        childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, shell=True)
        result, error_message = childProcess.communicate()
        if six.PY3:
            result = result.decode()
            error_message = error_message.decode()
        print(result.strip())
        try:
            if cluster_type == 'banting':
                jobnumber = int(result.strip().split('.')[0])
            elif cluster_type == 'beagle':
                jobnumber = int(result.strip().split(' ')[2])
            elif cluster_type is None:
                jobnumber = -1
        except Exception as e:
            print(e)
            print('Error! Retrying...')
        if n_tries > 5:
            time.sleep(60)
    return jobnumber


#%%
def get_running_jobs():
    """Get a list of job ids for jobs that are still running.
    """
    result = _get_running_jobs()
    running_jobs = _parse_running_jobs(result)
    return running_jobs


def _get_running_jobs():
    """Use qstat to get all running jobs for a user (in this case "strokach")
    """
    while True:
        child_process = subprocess.Popen('qstat | grep strokach',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, error_message = child_process.communicate()
        rc = child_process.returncode
        if six.PY3:
            result = result.decode()
            error_message = error_message.decode()
        if error_message and rc:
            print('An error occured while trying to get a list of running jobs!\n{};{};{};\n'
                  .format(result, error_message, rc))
            time.sleep(60)
            continue
        break
    return result


def _parse_running_jobs(result):
    """Parse the results produced by ``_get_running_jobs`` function.
    """
    running_jobs = []
    for line in result.split('\n'):
        try:
            jobnumber = int(line.strip().split()[0].split('.')[0])
            running_jobs.append(jobnumber)
        except (IndexError, ValueError):
            print(line)
            print('Strange line above, skipping...')
    return running_jobs



def _run_system_command_safely(system_command, max_num_of_tries=10):
    """
    """
    number_of_tries = 0
    while True:
        number_of_tries += 1
        child_process = subprocess.Popen(
            system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = result.decode()
            error_message = error_message.decode()
        # Rerun command if you got an error message and you tried < `max_num_of_tries` times
        if error_message and number_of_tries < max_num_of_tries:
            print(error_message)
            time.sleep(60)
        else:
            break
    if error_message:
        raise Exception('Quitting because of unresolvable errors: {}'.format(error_message))
    return result


def get_num_submitted_jobs():
    """Use this if you don't need job ids for every running job
    """
    system_command = 'jobstatus | grep strokach'
    result = _run_system_command_safely(system_command)
    jobstatus = result.strip().split('\n')
    num_submitted_jobs = sum([sum([int(x) for x in js.split()[2:]]) for js in jobstatus])
    return num_submitted_jobs


def get_num_submitted_jobs_by_all():
    """Returns the total number of jobs running on the cluster
    """
    system_command = 'jobstatus | grep totals'
    result = _run_system_command_safely(system_command)
    jobstatus = result.strip().split()
    if jobstatus[0] != 'totals':
        raise Exception('Unexpected `jobstatus` output: {}'.format(jobstatus))
    total = sum(int(x) for x in jobstatus[1:])
    return total


#%%
def wait_for_jobs(running_jobs):
    """Wait untill all the jobs in the `job_ids` list are done.
    """
    if cluster_type is None:
        return
    while running_jobs:
        time.sleep(60)
        running_jobs = set(get_running_jobs())


