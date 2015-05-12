from optparse import OptionParser;
import sys;
import zipfile;

parser = OptionParser();

def main(arguments):
    (options, args) = parser.parse_args(arguments);
    assert( len(args) > 0 );
    with zipfile.ZipFile(args[0], 'r') as zip:
		zip.extractall();
    return 0;
	
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]));

