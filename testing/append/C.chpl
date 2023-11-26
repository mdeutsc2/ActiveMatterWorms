module C {
    extern {
        #include <stdio.h>

        static int appendFileExists(const char *filename) {
            FILE *file = fopen(filename, "r");

            if (file != NULL) {
                fclose(file);
                return 1; // File exists
            }

            return 0; // File does not exist
        }

        static int appendToFile(const char *filename, const char *str) {
            FILE *file = fopen(filename, "a");

            if (file == NULL) {
                perror("Error opening file");
                return 1;
            }

            fprintf(file, "%s\n", str);

            fclose(file);
            return 0;
        }
    }
}