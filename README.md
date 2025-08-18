// created by shazia
// first mini-task of bioinformatics 1 in 4th semester which is now slightly updated by me
// it Reads two strings, builds the DP matrix, traces back the best local path, and prints the matrix, best score, and aligned substrings.
// Default scoring: match = +2, mismatch = −1, gap = −1


#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cstdint>

struct SWResult {
    int bestScore{};
    std::string aln1;
    std::string aln2;
    std::vector<std::vector<int>> H;
// score matrix
};

enum class Dir : std::uint8_t { Stop = 0, Diag = 1, Up = 2, Left = 3 };

// Smith–Waterman local alignment
SWResult smith_waterman(const std::string& a,
                        const std::string& b,
                        int match = 2,
                        int mismatch = -1,
                        int gap = -1)
{
    const int n = static_cast<int>(a.size());
    const int m = static_cast<int>(b.size());

    std::vector<std::vector<int>>  H(n + 1, std::vector<int>(m + 1, 0));
    std::vector<std::vector<Dir>> dir(n + 1, std::vector<Dir>(m + 1, Dir::Stop));

    int bestScore = 0, bestI = 0, bestJ = 0;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const int s = (a[i - 1] == b[j - 1]) ? match : mismatch;

            const int d = H[i - 1][j - 1] + s; // diag
            const int u = H[i - 1][j] + gap;   // up (gap in b)
            const int l = H[i][j - 1] + gap;   // left (gap in a)

            int cell = 0;
            Dir choice = Dir::Stop;

            if (d >= u && d >= l && d > 0)      { cell = d; choice = Dir::Diag; }
            else if (u >= d && u >= l && u > 0) { cell = u; choice = Dir::Up;   }
            else if (l >= d && l >= u && l > 0) { cell = l; choice = Dir::Left; }
            // else keep 0/Stop for local alignment

            H[i][j] = cell;
            dir[i][j] = choice;

            if (cell > bestScore) { bestScore = cell; bestI = i; bestJ = j; }
        }
    }

    // Traceback
    std::string aln1, aln2;
    int i = bestI, j = bestJ;
    while (i > 0 && j > 0 && dir[i][j] != Dir::Stop) {
        switch (dir[i][j]) {
            case Dir::Diag:
                aln1.push_back(a[i - 1]);
                aln2.push_back(b[j - 1]);
                --i; --j;
                break;
            case Dir::Up:
                aln1.push_back(a[i - 1]);
                aln2.push_back('-');
                --i;
                break;
            case Dir::Left:
                aln1.push_back('-');
                aln2.push_back(b[j - 1]);
                --j;
                break;
            default: 
            // safety
                i = 0; j = 0;
                break;
        }
    }
    std::reverse(aln1.begin(), aln1.end());
    std::reverse(aln2.begin(), aln2.end());

    return { bestScore, aln1, aln2, std::move(H) };
}

static void print_matrix(const std::vector<std::vector<int>>& H,
                         const std::string& a, const std::string& b)
{
    const int n = static_cast<int>(a.size());
    const int m = static_cast<int>(b.size());

    std::cout << "\nScore matrix (H):\n    ";
    for (int j = 0; j < m; ++j) std::cout << std::setw(3) << b[j];
    std::cout << "\n";

    for (int r = 0; r <= n; ++r) {
        if (r == 0) std::cout << "  ";
        else        std::cout << a[r - 1] << " ";
        for (int c = 0; c <= m; ++c) std::cout << std::setw(3) << H[r][c];
        std::cout << "\n";
    }
}

int main() {
    std::string base1, base2;
    std::cout << "Enter first string: ";
    std::cin  >> base1;
    std::cout << "Enter second string: ";
    std::cin  >> base2;

    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);


    constexpr int MATCH = 2;
    constexpr int MISMATCH = -1;
    constexpr int GAP = -1;

    SWResult r = smith_waterman(base1, base2, MATCH, MISMATCH, GAP);

    print_matrix(r.H, base1, base2);

    std::cout << "\nBest local alignment score: " << r.bestScore << "\n";
    std::cout << r.aln1 << "\n";
    for (std::size_t k = 0; k < r.aln1.size(); ++k) {
        if (r.aln1[k] == r.aln2[k]) std::cout << "|";
        else if (r.aln1[k] == '-' || r.aln2[k] == '-') std::cout << " ";
        else std::cout << ".";
    }
    std::cout << "\n" << r.aln2 << "\n";
    return 0;
}
