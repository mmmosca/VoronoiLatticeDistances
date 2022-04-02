#ifndef CMD_H_
#define CMD_H_

#ifdef DEBUG
#define PRINT_DATA()	printf("CURRENT ---- OPTIND: %d ---- ARGVIND: %d ---- OPTION: %c ---- OPTARG: %s\n", optind, argvind, argv[optind][1], myoptarg);
#define PRINT_DATAW()	printf("CURRENT ---- OPTIND: %d ---- ARGVIND: %d ---- OPTION: %s ---- OPTARG: %s\n", optind, argvind, argv[optind], myoptarg);
#else
#define PRINT_DATA()	printf("OPTION: %c ---- OPTARG: %s\n", argv[optind][1], myoptarg);
#define PRINT_DATAW()	printf("OPTION: %s ---- OPTARG: %s\n", argv[optind], myoptarg);
#endif

#include <stdio.h>
#include <string.h>
#include <Windows.h>

int isCharInString(char c, char* str);

int AreStringsEqualFrom(const char* s1, const char* s2, int from);

int isSubstring(char* sub, char* str);

char* strsep(char** elem_pointer, char* pattern);

struct CommandLine {
private:
	/*option index in the array of line commands*/
	/*argument index in the array of line commands*/
	int optind, argvind, start, formatind, argformatind;
	/*pointer to the argument*/
	char *curr_option;
public:
	char* myoptarg;

	CommandLine() : start{ 1 }, formatind{ -1 }, argformatind{ 0 } {
		reset_values();
	};

	void reset_values();

	char mygetopt(int argc, char** argv, char* format);

	char* mygetoptW(int argc, char** argv, char* format);
};
#endif /* CMD_H_ */