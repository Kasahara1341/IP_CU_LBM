#pragma once

class Point {
public:
    double x;
    double y;

    // コンストラクタ
    Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {}

    // メンバ関数
    double magnitudea() const {
        return std::sqrt(x*x + y*y);
    }
};
