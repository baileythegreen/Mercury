class EnsemblAPIError(ValueError):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

    def report(self, message = "Something's wrong with: "):
        print("%s%s" % (message, self.value))


class EnsemblClientError(EnsemblAPIError):
    def report(self):
        super().report("There's a problem setting up the client: ")

class EnsemblRequestError(EnsemblAPIError):
    def report(self):
        super().report("There's a problem with the request: ")

class EnsemblFileError(EnsemblAPIError):
    def report(self):
        super().report("There's a problem with this file: ")

## client object
class NonexistantGenomeBuildError(EnsemblClientError):
    def report(self):
        super().report("This build does not seem to exist: ")

class HeaderIsNotDictError(EnsemblClientError):
    def report(self):
        super().report("This header is not a dict: ")

class ExtensionIsNotStringError(EnsemblClientError):
    def report(self):
        super().report("This extension is not a string: ")

class EnsemblServerRequestError(EnsemblClientError):
    def report(self):
        super().report("Something has gone wrong with the request; check the url, headers, and extensions are valid: ")



## request object
class UnsupportedFeaturesContainerError(EnsemblRequestError):
    def report(self):
        super().report("This features container should be a list, set, or tuple: ")

class FeatureIsNotStringError(EnsemblRequestError):
    def report(self):
        super().report("This feature is not a string: ")

class SublevelTooManyPartsError(EnsemblRequestError):
    def report(self):
        super().report("This tuple has too many elements; expected a parent string and a child list, in that order: ")

class SublevelTooFewPartsError(EnsemblRequestError):
    def report(self):
        super().report("This tuple has too few elements; expected a parent string and a child list, in that order: ")

class SublevelNotInTupleFormError(EnsemblRequestError):
    def report(self):
        super().report("Sublevel should be a tuple containing a parent string and a child list, in that order ")

class SnpIdsNotListError(EnsemblRequestError):
    def report(self):
        super().report("SNP ids must be in list form to keep requests.post() from complaining: ")

class SomethingInRequestIsEscapedError(EnsemblRequestError):
    def report(self):
        super().report("There is something, probably a quote character escaped and that's a problem: ")

class DecodedRecordsAreNotListError(EnsemblRequestError):
    def report(self):
        super().report("The decoded records are not a list and tey should be: ")

class UnrecognisedFeatureNameError(EnsemblRequestError):
    def report(self):
        super().report("This feature name does not appear in the decoded records: ")

class RecordHasNoIdError(EnsemblRequestError):
    def report(self):
        super().report("This entry does not have an id field, which is weird: ")

class FieldsNotRequestedFeaturesObjectError(EnsemblRequestError):
    def report(self):
        super().report("A RequestedFeatures object is required as input; RequestedFeatures(toplevel, sublevel):")


## file reader
class FilenameNotStringError(EnsemblFileError):
    def report(self):
        super().report("Filename must be a string: ")

class JsonFileFormatError(EnsemblFileError):
    def report(self):
        super().report("It looks like the input JSON file is malformed: ")

class JsonIsEmptyError(EnsemblFileError):
    def report(self):
        super().report("It looks like the input JSON file is empty.")

class NoAnalysesFoundError(EnsemblFileError):
    def report(self):
        super().report("It looks like there are no analyses in the input JSON file.")

class NoSNPsToQueryError(EnsemblFileError):
    def report(self):
        super().report("It looks like none of the analyses had GWAS hits.")

class ThisIsNotAFilenameError(EnsemblFileError):
    def report(self):
        super().report("You can pass a filename here, but this isn't one: ")

# Wrong data types
# class IsNotStringError(EnsemblClientError):
    # def report(error):
        # print("This value is not a string: ")
