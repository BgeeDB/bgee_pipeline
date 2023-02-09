import json
import fileinput
import io
import shutil

from typing.io import TextIO


class Util:
    """ Utility class with general purpose methods."""

    @staticmethod
    def convert_result_to_key_value(column_key: str, results: object) -> dict:
        """Return the query tabular results as a key-value data structure where the key are the values of a given column
         (column_key).

        :param column_key: the column value that will be used as a key in the returned dictionary.
        :param results: the results of a SPARQL query as returned by the Object(SPARQLWrapper).query().convert() method.
        :return:
             a dictionary.
        """
        # get header (column names) from results
        header = results["results"]["bindings"][0].keys()
        # display table of results:
        table = {}
        # the SPARQL JSON results to the query are available in the "results", "bindings" entry:
        for entry in results["results"]["bindings"]:
            # append entries from the results to a key-value dictionary where the key is a cell value in column_key
            row = {}
            for column in header:
                if column != column_key:
                    row.update({entry[column_key]["value"]: entry[column]["value"]})
                else:
                    row.update({entry[column_key]["value"]: ""})
            table.update(row)
        return table

    @staticmethod
    def rewrite_results_csv(output_header: str, file_path: str, results: object,
                            columns_to_replace: list = None, mapping: dict = None,
                            drop_duplicates: bool = False):
        """Rewrite JSON-like SPARQL results of a projected variable (i.e. column_to_replace)
         based on a mapping dictionary. If no column_to_replace is provided then this method will write a CSV file
         with the results without performing any mapping.

        :param output_header: The string header of the output files where columns should be comma separated.
        :param file_path: The full file directory path.
        :param results: the results of a SPARQL query as returned by the Object(SPARQLWrapper).query().convert(JSON)
         method.
        :param columns_to_replace: columns where we want to replace all values according to a mapping dictionary.
            Default is None.
        :param mapping: a dictionary to map all cell values in the column_to_replace to a corresponding one
         defined in this mapping. Default is None.
        """
        header = results["results"]["bindings"][0].keys()
        output = io.StringIO()
        output.write(output_header + "\n")
        file = open(file_path, 'a')
        seen_set = set()
        newline_count = 1
        for entry in results["results"]["bindings"]:
            line = ""
            mapped_value = None
            for column in header:
                if columns_to_replace is not None and column in columns_to_replace:
                    try:
                       cell_value = int(entry[column]["value"])
                    except:
                       cell_value = str(entry[column]["value"])
                    mapped_value = mapping.get(cell_value)
                    if mapped_value is None:
                        line += entry[column]["value"] + ','
                    else:
                        if isinstance(mapped_value, list) and len(mapped_value) == 1:
                            line += str(mapped_value[0]) + ','
                        else:
                            line += str(mapped_value) + ','
                else:
                    line += entry[column]["value"] + ','
            line = line[:-1]
            if mapped_value is not None and isinstance(mapped_value, list) and len(mapped_value) > 1:
                line_list = []
                for value in mapped_value:
                    line_list.append(line.replace(str(mapped_value), str(value)))
            else:
                line_list = [line]
            for line in line_list:
                if drop_duplicates:
                    if line in seen_set:  # skip duplicate
                        continue
                    seen_set.add(line)
                output.write(line + '\n')
                if newline_count % 10000 == 0:
                    output.seek(0)
                    shutil.copyfileobj(output, file)
                    output.close()
                    output = io.StringIO()
                newline_count = newline_count + 1
        output.seek(0)
        shutil.copyfileobj(output, file)
        output.close()
        file.close()

    @staticmethod
    def replaces_item_string_list(string_list: list, old_item: object, new_item: object):
        """Replaces an item from a string list.

        :param string_list: a string list.
        :param old_item: The object with a string representation to be replaced.
        :param new_item: The replacement object with a string representation.
        """
        new_item = str(new_item)
        old_item = str(old_item)
        if old_item in string_list:
            string_list.remove(old_item)
            string_list.append(new_item)

    @staticmethod
    def drop_duplicates(file_path: str):
        """

        :param file_path:
        """
        seen = set()  # set for fast O(1) amortized lookup
        for line in fileinput.FileInput(file_path, inplace=1):
            if line in seen: continue  # skip duplicate
            seen.add(line)
            print(line.strip('\n'), )  # standard output is now redirected to the file

    @staticmethod
    def get_asymmetric_pairs_from_list(string_list: list, include_same_item_pairs: bool) -> list:
        """ Create the item pair list excluding symmetric pairs such as (a,b) = (b,a),
        thus only one of these pairs is kept.

        :param string_list: A list of strings
        :param include_same_item_pairs: If True, it includes pairs with the same species in the output list.
        :return:
            A list of binary-tuples, e.g. [(a,b),(a,c)].
        """
        values_species = []
        for entry in string_list:
            for entry2 in string_list:
                if entry != entry2 or include_same_item_pairs:
                    if not (entry2, entry) in values_species:
                        values_species.append((entry, entry2))
        return values_species

    @staticmethod
    def set_as_first_in_pair(pair_list: list, items: list) -> list:
        """ Set an item from a pair as first based on a list of "must be first" items. If all pair elements are in
        a given items list, there is no changes in the pair order.  

        :param pair_list: A list of binary-tuples, e.g. [(a,b),(a,c)]
        :param items: possible items to set as first in a pair. If all pair elements are in
        this list, there is no changes in the pair order.
        :return:
            A list of binary-tuples, e.g. [(a,b),(a,c)].
        """
        for entry1, entry2 in pair_list:
            if entry2 in items and entry1 not in items:
                pair_list.remove((entry1, entry2))
                pair_list.append((entry2, entry1))
        return pair_list


    @staticmethod
    def remove_pairs(pairs_to_remove_file_path: str, pairs: list):
        """Remove pairs from a list of pairs.

        :param pairs_to_remove_file_path: the path of a JSON file that contains the pairs to be removed.
            Example of this file contents: [(a,b),(a,c)]
        :param pairs: the list of pairs where a pair is a binary tuple. Example: [(a,b),(a,c),(b,c)].
        """
        try:
            with open(pairs_to_remove_file_path, "r") as file:
                pairs_remove_list = Util.read_tmp_file(file)
                try:
                    for [item1, item2] in pairs_remove_list:
                        if pairs.__contains__((item1, item2)):
                            pairs.remove((item1, item2))
                        else:
                            if pairs.__contains__((item2, item1)):
                                pairs.remove((item2, item1))
                except:
                    raise ValueError(
                        "Syntax error in " + pairs_to_remove_file_path + " file, it should be a list of lists"
                        + ", e.g. [[1,2],[3,4]].")
        except FileNotFoundError as file_error:
            raise FileNotFoundError("The file {} does not exist.".format(pairs_to_remove_file_path))

    @staticmethod
    def read_tmp_file(tmp_file: TextIO):
        file_content = tmp_file.read().__str__()
        try:
            pairs_remove_list = list(json.loads(file_content))
        except:
            try:
                pairs_remove_list = list(json.loads(
                    "[" + file_content.strip(',') + "]"))
            except json.JSONDecodeError as e:
                raise e
        return pairs_remove_list


    @staticmethod
    def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
        # Print New Line on Complete
        if iteration == total:
            print()

