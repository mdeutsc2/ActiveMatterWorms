#include <stdio.h>

static void appendFile(const char *filename) {
    FILE *file = fopen(filename, "a"); // Open the file in append mode

    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    for (int i = 1; i <= 100; ++i) {
        fprintf(file, "%d\n", i); // Append each number to the file
    }

    fclose(file); // Close the file
}

/*
int main() {
    const char *filename = "numbers.txt";
    appendNumbersToFile(filename);

    return 0;
}
*/