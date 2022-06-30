"""
Execute binary with argument. If the argument is pattern, then start it with $.
"""

from subprocess import Popen, PIPE
import sys
import os
import glob


def execute_binary(args):
    process = Popen(' '.join(args), shell=True, stdout=PIPE, stderr=PIPE)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode

    return rc, std_output, std_error

# generate execution arguments
def generate_args(binary, *params):
    arguments = [binary]
    arguments.extend(list(params))
    return arguments


def check_correctness(binary_path, data_graph_path, data_graph_update_path, query_folder_path, expected_results):
# find all query graphs.
    query_graph_path_list = glob.glob('{0}/*'.format(query_folder_path))

    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path,
                                       '-u', data_graph_update_path, '-num', 'MAX')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if 'Result count' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break

            expected_embedding_num = expected_results[query_graph_name]
            if embedding_num != expected_embedding_num:
                print('Query {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name, expected_embedding_num,
                                                                  embedding_num))
                exit(-1)

            print('pass {0}:{1}'.format(query_graph_name, embedding_num))
        else:
            print('Query {0} execution error.'.format(query_graph_name))
            exit(-1)


if __name__ == '__main__':
    input_binary_path = sys.argv[1]
    if not os.path.isfile(input_binary_path):
        print('The binary {0} does not exist.'.format(input_binary_path))
        exit(-1)

    # load expected results.

    dir_path = os.path.dirname(os.path.realpath(__file__))
    input_data_graph_path = '{0}/insert/data_graph/data.graph'.format(dir_path)
    input_data_graph_update_path = '{0}/insert/data_graph/insertion.graph'.format(dir_path)
    input_query_graph_folder_path = '{0}/insert/query_graph/'.format(dir_path)
    input_expected_results_file = '{0}/insert/expected_output.txt'.format(dir_path)

    input_expected_results = {}

    with open(input_expected_results_file, 'r') as f:
        for line in f:
            if line:
                result_item = line.split(':')
                input_expected_results[result_item[0].strip()] = int(result_item[1].strip())

    check_correctness(input_binary_path, input_data_graph_path, input_data_graph_update_path,
                      input_query_graph_folder_path, input_expected_results)

    print('Pass all insert test cases.')

    input_data_graph_path = '{0}/delete/data_graph/data.graph'.format(dir_path)
    input_data_graph_update_path = '{0}/delete/data_graph/deletion.graph'.format(dir_path)
    input_query_graph_folder_path = '{0}/delete/query_graph/'.format(dir_path)
    input_expected_results_file = '{0}/delete/expected_output.txt'.format(dir_path)

    input_expected_results = {}

    with open(input_expected_results_file, 'r') as f:
        for line in f:
            if line:
                result_item = line.split(':')
                input_expected_results[result_item[0].strip()] = int(result_item[1].strip())

    check_correctness(input_binary_path, input_data_graph_path, input_data_graph_update_path,
                      input_query_graph_folder_path, input_expected_results)
    print('Pass all delete test cases')