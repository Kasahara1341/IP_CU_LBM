#include "all.hpp"

// chatGPT製
// fileから指定列(column_index)を読み込んでvecに格納する
// column_index は 0始まり (1列目を使いたければ0, 2列目なら1)
void readCSV(const std::string& filename, std::vector<double>& vec, int column_index) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string token;
        int col = 0;

        bool success = false;  // 数値に変換できたかどうか
        double value = 0.0;

        while (std::getline(ss, token, ',')) {
            if (col == column_index) {
                try {
                    value = std::stod(token);  // 数値に変換できるか試す
                    success = true;
                } catch (const std::invalid_argument&) {
                    success = false;  // 数値変換失敗 → ヘッダーや文字列行
                }
                break;  // 目的の列だけ取得したらループ終了
            }
            col++;
        }

        if (success) {
            vec.push_back(value);  // 数値に変換できたら格納
        }
        // 数値にならなければ次の行へ
    }
}
// chatGPT製
// 2次元ベクトル(行：時間，列：場所)にまとめたデータをCSVに保存する関数
void writeCSV(const std::string &filename,
              const std::vector<std::vector<double>> &data) {
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "ファイルを開けませんでした: " << filename << "\n";
        return;
    }

    for (size_t t = 0; t < data.size(); t++) {
        for (size_t i = 0; i < data[t].size(); i++) {
            ofs << data[t][i];
            if (i != data[t].size() - 1) ofs << ",";
        }
        ofs << "\n";
    }

    ofs.close();
    std::cout << "CSV出力完了: " << filename << "\n";
}
