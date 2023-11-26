use IO;
use IO.FormattedIO;
use CTypes;
import C;
//require "appendFile.h", "appendFile.c";

// prototype caller for external C function
//extern proc appendFile(str:c_ptrConst(c_char));

proc main() {
    var filename:string = "numbers.txt";
    writeln("testing append C-extension");
    var ret_int:int, test_str:string;
    if C.appendFileExists(filename.c_str()) != 0 {
        var error_str:string = filename+" already exists";
        halt(error_str);
    }
    for i in 1..5 {
        test_str = i:string;
        ret_int = C.appendToFile(filename.c_str(),test_str.c_str());
        if ret_int > 0 {
            halt("error in appending");
        }
    }
}