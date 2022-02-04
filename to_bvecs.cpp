#include "cmdline.h"
#include "misc.hpp"

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);

    cmdline::parser p;
    p.add<std::string>("input_fn", 'i', "input file name of database sketches in ascii format", true);
    p.add<std::string>("output_fn", 'o', "output file name of database sketches in bvec format", true);
    p.parse_check(argc, argv);

    auto input_fn = p.get<std::string>("input_fn");
    auto output_fn = p.get<std::string>("output_fn");

    std::ifstream ifs(input_fn);
    if (!ifs) {
        std::cerr << "open error: " << input_fn << '\n';
        return 1;
    }

    std::ofstream ofs(output_fn);
    if (!ofs) {
        std::cerr << "open error: " << output_fn << '\n';
        return 1;
    }

    std::vector<uint8_t> vec;
    for (std::string line; std::getline(ifs, line);) {
        vec.clear();
        std::istringstream iss(line);
        for (std::string s; iss >> s;) {
            auto v = std::stoul(s);
            if (256 <= v) {
                std::cerr << "error: input value must be < 256: " << v << '\n';
                return 1;
            }
            vec.emplace_back(static_cast<uint8_t>(v));
        }
        uint32_t dim = vec.size();
        ofs.write(reinterpret_cast<const char*>(&dim), sizeof(uint32_t));
        ofs.write(reinterpret_cast<const char*>(vec.data()), dim);
    }

    ifs.close();
    ofs.close();

    return 0;
}
