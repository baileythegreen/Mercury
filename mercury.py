import requests, sys, os
import pprint
import json
import optparse
from EnsemblAPIErrors import *
import ast


def doublequoter(s):
    return json.dumps(ast.literal_eval(s))

class EnsemblClient():
    def __init__(self, build = None, headers = None, extension = None):
        self.server = build
        self.headers = headers
        self.extension = extension

    @property
    def server(self):
        return self.__server

    @server.setter
    def server(self, value):
        known_builds = {"37" : "grch37.", "38" : ""}
        try:
            if value not in known_builds:
                raise NonexistantGenomeBuildError(value)
            else:
                self.__server = "https://%srest.ensembl.org" % known_builds[value]
        except NonexistantGenomeBuildError as error:
            error.report()

    @property
    def headers(self):
        return self.__headers

    @headers.setter
    def headers(self, value):
        try:
            if not isinstance(value, dict):
                raise HeaderIsNotDictError(value)
            else:
                self.__headers = value
        except HeaderIsNotDictError as error:
            error.report()

    @property
    def extension(self):
        return self.__extension

    @extension.setter
    def extension(self, value):
        try:
            if not isinstance(value, str):
                raise ExtensionIsNotStringError(value)
            else:
                self.__extension = value
        except ExtensionIsNotStringError as error:
            error.report()

    @property
    def snp_list(self):
        return self.__snp_list

    @snp_list.setter
    def snp_list(self, value):
        try:
            if not isinstance(value, list):
                raise SnpIdsNotListError(value)
            else:
                self.__snp_list = doublequoter("{ 'ids' : %s }" % value)
        except SnpIdsNotListError as error:
            error.report()

### non-setter methods
    def make_request(self, value = None):
        ## value = self.snp_list_request
        try:
            if '\\' in self.snp_list:
                raise SomethingInRequestIsEscapedError(value)
            records = requests.post(self.server + self.extension,
                headers = self.headers, data = value)
            if not records.ok:
                records = requests.post(self.server + self.extension,
                    headers = self.headers, json = value)
            if not records.ok:
                raise EnsemblServerRequestError(records.raise_for_status())
        except EnsemblServerRequestError as error:
            error.report()
        return records

    def decode_return(self, records, outfile = None):
        # this should write to a file, or something, but the file name should be unique?
        decoded = records.json()
        return decoded

class Extractor():
    def __init__(self, data = None, fields = None):
        self.data = data  # file-like object, or input stream
        self.fields = fields

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data):
        try:
            if data is None:
                self.__data = data
            elif isinstance(data, str) and os.path.isfile(data):
                __data = ast.literal_eval(open(data).read())
                if not isinstance(__data, list):
                    raise ThisIsNotAFilenameError(data)
                else:
                    data = __data
            elif not isinstance(data, list):
                raise DecodedRecordsAreNotListError(data)
            else:
                self.__data = data
        except EnsemblRequestError as error:
            error.report()

    @property
    def fields(self):
        return self.__fields

    @fields.setter
    def fields(self, fields):
        try:
            if fields is None:
                self.__fields = fields
            elif not isinstance(fields, RequestedFeatures):
                raise FieldsNotRequestedFeaturesObjectError(fields)
            else:
                self.__fields = fields
        except EnsemblRequestError as error:
            error.report()

    def extract_info(self, data = None, fields = None):
        def extract_snp_results(snp):
            result = {}
            for feature in self.fields.toplevel:
                try:
                    if isinstance(feature, str):
                        result[feature] = snp[feature]
                    elif isinstance(feature, tuple):
                        parent, children = feature
                        result[parent] = {child: snp[parent][0][child] for child in children}
                except KeyError: # this should not prevent the code from continuing
                    self.not_found(feature)
                    continue
            return result
        try:
            if isinstance(self.data, list):
                return [extract_snp_results(record) for record in self.data]
            else:
                raise DecodedRecordsAreNotListError(self.data)
        except EnsemblRequestError as error:
            error.report()

        # definition of behaviour triggered by an unfound feature
        # could be a coroutine
    def not_found(self, feature):
        pass

class RequestedFeatures():
    def __init__(self, toplevel = None, sublevel = None):
        self.sublevel = self.set_sublevel(sublevel)
        self.toplevel = self.set_toplevel(toplevel)

    def __init__(self, toplevel = None, sublevel = None):
        self.sublevel = self.set_sublevel(sublevel)
        self.toplevel = self.set_toplevel(toplevel)

    def set_toplevel(self, value):
        try:
            if isinstance(value, str):
                feats = [value]
            elif isinstance(value, (list, set, tuple)):
                feats = []
                for feature in value:
                    try:
                        if isinstance(feature, tuple):
                            subfeatures = self.set_sublevel(feature)
                            feats.append(subfeatures)
                        elif isinstance(feature, str):
                            feats.append(feature)
                        else:
                            raise FeatureIsNotStringError(feature)
                    except FeatureIsNotStringError as error:
                        print('This feature is not a string: %s' % error)
            else:
                raise UnsupportedFeaturesContainerError(value)
        except EnsemblRequestError as error:
            error.report()
        return feats

    def set_sublevel(self, value):
        try:
            if isinstance(value, tuple):
                num_parts = len(value)
                if num_parts == 2:
                    subfeats = []
                    parent, children = value
                    if not isinstance(parent, str):
                        raise FeatureIsNotStringError(parent)
                    for child in children:
                        if isinstance(child, str):
                            subfeats.append(child)
                        else:
                            raise FeatureIsNotStringError(child)
                else:
                    if len(value) > 2: raise SublevelTooManyPartsError(value) ###### an error saying it's too long
                    elif len(value) < 2: raise SublevelTooFewPartsError(value)
                    else: raise EnsemblRequestError(tuple)
                return (parent, subfeats)
            else:
                raise SublevelNotInTupleFormError(value)
        except EnsemblRequestError as error:
            error.report()

class InputFileReader():
    def __init__(self, filename = None):
        self.file_toread = filename if self.check_filename(filename) else None
        self.loaded_data = self.read_file(self.file_toread)
        self.parsed_data = self.parse_dataset()

    def check_filename(self, filename):
        try:
            if not isinstance(filename, str): raise FilenameNotStringError(filename)
            else: return True
        except EnsemblAPIError as error:
            error.report()
            return False

    def read_file(self, filename):
        try:
            with open(filename, 'r') as f:
                data = json.load(f)
        except EnsemblAPIError as error:
            error.report(error)
        return data

    def parse_dataset(self):
        parsed_data = {}
        subcohorts = []
        snp_list = []
        analyses_list = []
        try:
            while self.loaded_data:
                subcohorts.append(self.loaded_data.popitem())
            if not subcohorts:
                raise JsonIsEmptyError
            while subcohorts:
                group = subcohorts.pop()[1]
                while group:
                    k = group.pop()
                    analyses_list.append(k)
                    [snp_list.append(snp) for item in k.values() for snp in item]
            parsed_data.update({'subcohorts': subcohorts, 'analyses_list': analyses_list,
                        'snp_list': snp_list})
            if not analyses_list:
                raise NoAnalysesFoundError
            if not snp_list:
                raise NoSNPsToQueryError
        except EnsemblFileError as error:
            error.report()
        return parsed_data

class Printer():
    def __init__(self, outfile = None):
        return

    @property
    def outfile(self):
        return self.__outfile

    @outfile.setter
    def outfile(self, value):
        try:
            if not isinstance(value, str):
                raise FilenameNotStringError(value)
            else:
                self.__outfile = value
        except EnsemblFileError as error:
            error.report()

    def write(self, output):
        def write_out(entry, file):
            try:
                print(entry['id'], file)
                pprint.pprint(entry, stream = file)
            except KeyError:
                raise RecordHasNoIdError(entry)
        try:
            with open(self.outfile, 'a') as f:
                print
                [write_out(entry, f) for entry in output]
        except RecordHasNoIdError as error:
            error.report()

## possibly unnecessary
class SnpSet():
    def __init__(self):
        pass

    def set_snps(self):
        ### some validation code
        pass


#================================================================================

def help():
    print('''
    =======================================================================================================================
    Usage:

    python get_efo_disease_types.py -t <tree_file> -i <list_file>
    [, -o <out_file> [, -l <log_file> [, -s <supracategory>
    [, -c <categories> [, -d <duplicates> ] ] ] ] ]

    A tree file and a disease list file must be provided.

    An output file name and a log file name may be optionally specified;
    if these are not given, default values will be used:

    Output file default: ontoquery_results.json

    Log file default: ontoquery.log
 --------------------------------------------------------------------------
    For more information, run:
    python get_efo_disease_types.py -h
    or
    python get_efo_disease_types.py --help
    =======================================================================================================================
    ''')
    sys.exit(1)

def parse_options(values):
    """p is the parser that was defined in main()"""
    usage = "usage: %prog -i list_file [options]"
    version = "%prog 1.0"
    description = "Retrieves information from the Ensembl database"
    globals().update({'p' : optparse.OptionParser(usage = usage,
                      version = version, description = description)}) #,
                      #formatter = PrettyHelpFormatter())})
    p.add_option("-i", "--input", action = "store", dest = "list_file",
                help = "a .txt file containing a list of SNPs")
    p.add_option("-b", "--build", action = "store", dest = "build",
                help = "the genome build to use")
    p.add_option("-e", "--extension", action = "store", dest = "extension",
                help = "the extension to use")
    p.add_option("--headers", action = "store", dest = "headers",
                help = "the API headers to use")
    p.set_defaults(build = "37", extension = "/variation/homo_sapiens",
                headers = { "Content-Type" : "application/json", "Accept" : "application/json"},
                log_file = 'log.txt', work_dir = os.getcwd().lstrip())
    return p.parse_args()

##### main

def main():
    #if len(sys.argv) == 1: help()
    opts, args = parse_options(sys.argv)#'''
    print("opts: %s" % repr(opts))
    print("args: %s" % repr(args))
    trait, round = args[0].split("_")[0:2]
    snps = []
    with open(args[0], 'r') as f:
        lines = f.readlines()
        snps.extend([line.strip("\n") for line in lines])
        f.close()
    print(snps)
    # create client
    #request_fields = RequestedFeatures(toplevel = toplevel, sublevel = sublevel)
    clio = EnsemblClient(build = opts.build, headers = opts.headers, extension = opts.extension)
    #reader = InputFileReader(data_file)
    #clio.snp_list = reader.parsed_data['snp_list']
    clio.snp_list = snps
    records = clio.make_request(clio.snp_list)
    decoded = clio.decode_return(records)
    clio.extension = "/overlap/region/human/"
    extract = [(k, v['mappings'][0]['location'], v['MAF'], v['most_severe_consequence']) for k, v in decoded.items()]
    with open('%s_%s_SNP_summary.txt' % (trait, round), 'w') as f:
        for snp in extract:
            print('%s\t%s\t%s\t%s' % snp, file = f)
        f.close()
    locations = [snp[1] for snp in extract]
    loc_headers = [("content_type", "application/json"), ("accept", "application/json"), ("feature", "variation"), ("feature", "gene"), ("feature", "regulatory")]
    tail = ';'.join([k + '=' + v for k,v in loc_headers])
    results = []
    discard = ['- C', '- A', '- G', '- T', '', 'alleles: ', 'clinical_significance: []', 'clinical_significance: ', '- benign', '- risk factor']
    with open('%s_%s_overlap_results.txt' % (trait, round), 'w') as f:
        for location in locations:
            records = requests.get(clio.server + clio.extension + location + '?' + tail)
            if not records.ok:
                print("Could not retrieve data on: %s" % location, file = f)
                print("Attempted url: %s" % clio.server + clio.extension + location + '?' + tail, file = f)
                print('=' * 80, file = f)
                continue
            lines = records.text.split('\n')
            lines = [line.lstrip() for line in lines if line.startswith(' ') and line.lstrip() not in discard and not line.startswith('-')]
            results.append(lines)
        import re
        reformatted = []
        for item in results:
            for entry in item:
                print(entry, file = f)
            #d = []
            #for entry in item:
            #    print(entry)
            #    try:
            #        tup = tuple(entry.split(': '))
            #        d.append(tup)
            #    except:
            #        continue
            #reformatted.append(d)
            #print(repr(reformatted))
            #pprint.pprint(reformatted)
            #pprint.pprint([line.lstrip() for line in reformatted if line.startswith(' ') and line.lstrip() not in discard and not line.startswith('-')], stream = f)
            print('=' * 80, file = f)
        f.close()


if __name__ == '__main__':
    main()
