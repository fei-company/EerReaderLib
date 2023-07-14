#include <cstdlib>
#include <iostream>

#include "../ElectronCountedFramesDecompressor.h"

int main() {
	auto exampleName = "FoilHole_12572127_Data_8400797_8400799_20230213_134259_EER.eer";
	ElectronCountedFramesDecompressor decompressor(exampleName);

	std::cout << "Decompressor initialized!" << std::endl;

	return EXIT_SUCCESS;
}