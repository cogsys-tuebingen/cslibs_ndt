#include <ndt/grid/mask.hpp>

int main(int argc, char *argv[])
{
    ndt::grid::Mask<2> mask_bin_2;
    for(std::size_t i = 0 ; i < mask_bin_2.rows ; ++i) {
        for(std::size_t j = 0; j < mask_bin_2.cols ; ++j) {
            std::cout << mask_bin_2[i * mask_bin_2.cols + j] << " ";
        }
        std::cout << std::endl;
    }

    ndt::grid::Mask<3> mask_bin_3;
    for(std::size_t i = 0 ; i < mask_bin_3.rows ; ++i) {
        for(std::size_t j = 0; j < mask_bin_3.cols ; ++j) {
            std::cout << mask_bin_3[i * mask_bin_3.cols + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
