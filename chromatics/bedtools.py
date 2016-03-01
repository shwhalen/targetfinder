import io
import os
import pandas as pd
import subprocess
import tempfile

def read_bed(x, **kwargs):
    return pd.read_csv(x, sep = r'\s+', header = None, index_col = False, **kwargs)

def write_bed(df, fn, **kwargs):
    df = df.copy()
    df.iloc[:, [1, 2]] = df.iloc[:, [1, 2]].astype(int) # ensure coordinates are sorted as integer
    df.sort_values(df.columns.tolist()[:3], inplace = True)
    df.to_csv(fn, sep = '\t', header = False, index = False, **kwargs)

def sort_bed(x):
    if isinstance(x, pd.DataFrame):
        return x.sort_values(x.columns.tolist()[:3])
    raise Exception('Not a pandas DataFrame')

def bedtools(operation, left_input, right_input = None, left_names = None, right_names = None):
    # if first input is a dataframe, feed via stdin
    if isinstance(left_input, pd.DataFrame):
        left_input_fn = 'stdin'
    elif isinstance(left_input, str):
        left_input_fn = os.path.abspath(left_input)
    else:
        raise Exception('First input must be DataFrame or filename.')

    # if second input is a dataframe, write to a temp file to be removed later
    right_input_cleanup = False
    if isinstance(right_input, pd.DataFrame):
        right_input_fd, right_input_fn = tempfile.mkstemp()
        right_input_cleanup = True
        write_bed(right_input, right_input_fn)
    elif isinstance(right_input, str):
        right_input_fn = os.path.abspath(right_input)

    # create command line for one or two argument operations
    if operation.startswith('merge'):
        cmdline = 'bedtools {} "{}"'.format(operation, left_input_fn)
    else:
        cmdline = 'bedtools {} -a "{}" -b "{}"'.format(operation, left_input_fn, right_input_fn)

    # call bedtools, need shell = True for read permissions
    p = subprocess.Popen(cmdline, shell = True, stdout = subprocess.PIPE, stdin = subprocess.PIPE, close_fds = True)

    # if needed, write first dataframe in bed format to stdin
    if left_input_fn == 'stdin':
        left_input_buffer = io.StringIO()
        left_input.to_csv(left_input_buffer, sep = '\t', header = False, index = False)
        stdout, _ = p.communicate(input = left_input_buffer.getvalue().encode('utf-8'))
    else:
        stdout, _ = p.communicate()
    assert p.returncode == 0

    # if second input was a dataframe written to a temporary file, clean it up
    if right_input_cleanup:
        os.close(right_input_fd)
        os.remove(right_input_fn)
        assert not os.path.exists(right_input_fn)

    # infer column names
    if isinstance(left_input, pd.DataFrame) and left_names is None:
        left_names = left_input.columns.tolist()
    if isinstance(right_input, pd.DataFrame) and right_names is None:
        right_names = right_input.columns.tolist()

    if operation.startswith('intersect'):
        names = []
        if operation.find('-wa') != -1:
            names += left_names
        if operation.find('-wb') != -1:
            names += right_names
        if operation.find('-wo') != -1 or operation.find('-wao') != -1:
            names = left_names + right_names + ['overlap']
        if operation.find('-c') != -1:
            names = left_names + ['count']
        if operation.find('-u') != -1:
            names = left_names
        if operation.find('-loj') != -1:
            names = left_names + right_names
    elif operation.startswith('merge'):
        names = left_names
    elif operation.startswith('closest'):
        names = left_names + right_names
        if operation.find('-d') != -1:
            names.append('distance')
    else:
        names = left_names + right_names

    # create dataframe from bedtools output stored in stdout
    if len(stdout) == 0:
        return pd.DataFrame(columns = names)
    return read_bed(io.StringIO(stdout.decode('utf-8')), names = names)
