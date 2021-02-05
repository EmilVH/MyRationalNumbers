//
// Created by EmilVH on 27.10.2019.
//
#pragma once
#ifndef RATIONAL_H
#define RATIONAL_H

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <iostream>
#include <math.h>

typedef std::complex<long double> cld;
const long double PI = acosl(-1.0);

size_t findUpperDegreeOfTwo(size_t v) {
    size_t n = 1;
    while (n < v)
        n *= 2;
    return n;
}

size_t reverseBits(size_t i, size_t n) {
    size_t ans = 0;
    for (; n > 1; n /= 2) {
        ans = ans * 2 + (i % 2);
        i /= 2;
    }
    return ans;
}

typedef std::vector<cld>::iterator TIt;

void makeGeneralFFT(TIt begin, TIt end, cld q) {
    size_t n = end - begin;
    if (n == 1)
        return;

    for (size_t i = 0; i < n; ++i) {
        size_t revi = reverseBits(i, n);
        if (i < revi)
            std::swap(*(begin + i), *(begin + revi));
    }

    for (size_t l = 2; l <= n; l *= 2) {
        cld ql = q;
        for (size_t ll = n; ll > l; ll /= 2)
            ql *= ql;


        for (auto itbegin = begin; itbegin != end; itbegin += l) {
            cld w(1.0, 0.0);
            auto lit = itbegin;
            auto rit = itbegin + l / 2;
            auto eit = itbegin + l;
            for (; rit < eit; ++lit, ++rit) {
                cld u = *lit;
                cld v = w * (*rit);
                *lit = u + v;
                *rit = u - v;
                w *= ql;
            }
        }
    }
}

std::vector<cld> makeGeneralFFT(std::vector<cld> a, cld q) {
    makeGeneralFFT(a.begin(), a.end(), q);
    return a;
}


std::vector<cld> makeFFT(const std::vector<cld> &a) {
    long double ang = 2.0 * PI / a.size();
    return makeGeneralFFT(a, cld(cosl(ang), sinl(ang)));
}

std::vector<cld> makeInverseFFT(std::vector<cld> a) {
    long double ang = 2.0 * PI / a.size();
    a = makeGeneralFFT(a, cld(cosl(ang), -sinl(ang)));
    for (size_t i = 0; i < a.size(); ++i)
        a[i] /= a.size();
    return a;
}

template<typename T>
cld makeComplexNumber(T x) {
    return cld(x);
}

template<>
cld makeComplexNumber(long long x) {
    return cld(x, 0.0);
}

template<typename T>
std::vector<cld> makeComplexVector(const std::vector<T> &a, size_t n) {
    std::vector<cld> ca(n, cld(0.0, 0.0));
    for (size_t i = 0; i < a.size(); ++i)
        ca[i] = makeComplexNumber(a[i]);
    return ca;
}


template<typename T>
std::vector<T> makeTVector(const std::vector<cld> &a) {
    std::vector<T> ans(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        ans[i] = T(a[i]);
    return ans;
}

template<>
std::vector<long long> makeTVector(const std::vector<cld> &a) {
    std::vector<long long> ans(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        ans[i] = static_cast<long long>(floorl(a[i].real() + 0.5));
    return ans;
}

template<typename T>
std::vector<T> multiplicatePolynoms(const std::vector<T> &a,
                                    const std::vector<T> &b) {
    size_t n = findUpperDegreeOfTwo(std::max(a.size(), b.size()));
    n *= 2;

    std::vector<cld> ca = makeComplexVector(a, n);
    std::vector<cld> cb = makeComplexVector(b, n);

    ca = makeFFT(ca);
    cb = makeFFT(cb);

    for (size_t i = 0; i < n; ++i)
        ca[i] *= cb[i];

    std::vector<cld> cc = makeInverseFFT(ca);
    return makeTVector<T>(cc);
}


template<class curr>
void my_swap(curr &a, curr &b) {
    curr tmp = a;
    a = b;
    b = tmp;
}

template<class curr>
void my_reverse(curr first, curr last) {
    while ((first != last) && (first != --last)) {
        my_swap(*first++, *last);
    }
}

int min(const int a, const int b) {
    return a < b ? a : b;
}

struct BigInteger {
public:
    BigInteger();

    BigInteger(int a);

    BigInteger(const std::string &s);

    explicit operator bool() const;

    explicit operator int() const;

    std::string toString() const;

    int compare(const BigInteger &a, const BigInteger &b, bool byAbs) const;

    BigInteger abs() const;

    BigInteger &operator+=(const BigInteger &a);

    BigInteger &operator-=(const BigInteger &a);

    BigInteger &operator*=(const BigInteger &a);

    BigInteger &operator/=(const BigInteger &b);

    BigInteger &operator%=(const BigInteger &b);

    friend bool operator>(const BigInteger &a, const BigInteger &b);

    friend bool operator<(const BigInteger &a, const BigInteger &b);

    friend bool operator==(const BigInteger &a, const BigInteger &b);

    friend bool operator!=(const BigInteger &a, const BigInteger &b);

    friend bool operator<=(const BigInteger &a, const BigInteger &b);

    friend bool operator>=(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator+(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator-(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator*(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator/(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator%(const BigInteger &a, const BigInteger &b);

    BigInteger &operator++();

    BigInteger &operator--();

    BigInteger operator++(int);

    BigInteger operator--(int);

    BigInteger operator-() const;


private:
    std::vector<int> numb_;

    bool sign_ = true;

    static const int BASE_ = static_cast<int>(10);

    void deleteLeadingZeroes_();

    void addAbs_(BigInteger &a, const BigInteger &b);

    void substractAbs_(const BigInteger &first, const BigInteger &second);

    void div_(const BigInteger &a, const BigInteger &b, BigInteger &res, BigInteger &mod);

    int getRes_(int a, int b) const;

    void normalize_();
};

BigInteger::BigInteger() {
    numb_.push_back(0);
}

BigInteger::BigInteger(int a) {
    if (a < 0) {
        sign_ = false;
    }
    a = std::abs(a);
    if (a == 0) {
        numb_.push_back(0);
    } else {
        while (a > 0) {
            numb_.push_back(a % 10);
            a /= 10;
        }
    }
}

BigInteger::BigInteger(const std::string &s) {
    int strt = 0;
    if (s[0] == '-') {
        sign_ = false;
        strt++;
    }
    if (s == "0") {
        numb_.push_back(0);
    } else {
        for (int i = static_cast<int>(s.size() - 1); i >= strt; --i) {
            numb_.push_back(s[i] - '0');
        }
    }
}

void BigInteger::deleteLeadingZeroes_() {
    for (int i = static_cast<int>(numb_.size()) - 1; i >= 0; --i) {
        if (numb_[i] == 0) {
            numb_.pop_back();
        } else {
            break;
        }
    }
    if (numb_.empty()) {
        numb_.push_back(0);
    }
}

void BigInteger::addAbs_(BigInteger &a, const BigInteger &b) {
    int over = 0;
    size_t max_size = std::max(a.numb_.size(), b.numb_.size());
    for (size_t i = 0; i < max_size; ++i) {
        int first = 0;
        int second = 0;
        if (a.numb_.size() > i) first = a.numb_[i];
        if (b.numb_.size() > i) second = b.numb_[i];
        if (a.numb_.size() > i) {
            a.numb_[i] = first + second + over;
        } else {
            a.numb_.push_back(first + second + over);
        }
        over = 0;
        if (a.numb_[i] >= BASE_) {
            a.numb_[i] -= BASE_;
            over = 1;
        }
    }
    if (over > 0) {
        a.numb_.push_back(1);
    }

}

void BigInteger::substractAbs_(const BigInteger &first, const BigInteger &second) {
    int owed = 0;
    numb_.resize(std::max(first.numb_.size(), second.numb_.size()));
    for (size_t i = 0; i < second.numb_.size(); ++i) {
        if (first.numb_[i] < (second.numb_[i] + owed)) {
            numb_[i] = first.numb_[i] + BASE_ - (second.numb_[i] + owed);
            owed = 1;
        } else {
            numb_[i] = first.numb_[i] - (second.numb_[i] + owed);
            owed = 0;
        }
    }
    for (size_t i = second.numb_.size(); i < first.numb_.size(); ++i) {
        if (first.numb_[i] < owed) {
            numb_[i] = first.numb_[i] + BASE_ - (owed);
            owed = 1;
        } else {
            numb_[i] = first.numb_[i] - (owed);
            owed = 0;
        }
    }
    deleteLeadingZeroes_();

}

void BigInteger::div_(const BigInteger &a, const BigInteger &b, BigInteger &res, BigInteger &mod) {
    for (int i = static_cast<int>(a.numb_.size() - 1); i > -1; --i) {
        mod *= static_cast<int>(BASE_);
        mod += a.numb_[i];
        unsigned int l = 0;
        unsigned int r = BASE_;
        while (r - l > 1) {
            unsigned int m = (r + l) / 2;
            BigInteger tmp = b.abs() * static_cast<int>(m);
            if (tmp > mod) {
                r = m;
            } else {
                l = m;
            }
        }
        res.numb_.push_back(l);
        mod -= b.abs() * static_cast<int>(l);
    }
    my_reverse(res.numb_.begin(), res.numb_.end());
    res.sign_ = !(b.sign_ ^ a.sign_);
    mod.sign_ = !(b.sign_ ^ a.sign_);
    res.deleteLeadingZeroes_();
    mod.deleteLeadingZeroes_();
}

BigInteger::operator bool() const {
    return (*this) != 0;
}

BigInteger::operator int() const {
    int ans = 0;
    for (int i = static_cast<int>(numb_.size() - 1); i >= 0; --i) {
        ans *= BASE_;
        ans += numb_[i];
    }
    ans *= sign_;
    return ans;
}

std::string BigInteger::toString() const {
    std::string ans;
    if (!sign_ && compare(*this, 0, true) != 0) {
        ans.push_back('-');
    }
    for (size_t i = 0; i < numb_.size(); ++i) {
        ans.push_back(numb_[numb_.size() - 1 - static_cast<int>(i)] + '0');
    }
    return ans;
}

int BigInteger::compare(const BigInteger &a, const BigInteger &b, bool byAbs) const {
    if (byAbs) {
        if (a.numb_.size() != b.numb_.size()) {
            return getRes_(a.numb_.size(), b.numb_.size());
        }
        for (int i = static_cast<int>(a.numb_.size()) - 1; i >= 0; --i) {
            if (a.numb_[i] != b.numb_[i]) {
                return getRes_(a.numb_[i], b.numb_[i]);
            }
        }
        return 0;
    } else {
        if (a.sign_ == b.sign_) {
            return a.sign_ ? compare(a, b, true) : -compare(a, b, true);
        } else {
            return a.sign_ ? 1 : -1;
        }
    }
}

BigInteger BigInteger::abs() const {
    BigInteger res = *this;
    res.sign_ = true;
    return res;
}

BigInteger &BigInteger::operator+=(const BigInteger &a) {
    if (sign_ == a.sign_) {
        addAbs_(*this, a);
    } else {

        if (compare(*this, a, true) > -1) {
            substractAbs_(*this, a);
        } else {
            substractAbs_(a, *this);
            sign_ = !sign_;
        }
    }
    normalize_();
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &a) {

    if (sign_ == a.sign_) {
        if (compare(*this, a, true) > -1) {
            substractAbs_(*this, a);
        } else {
            substractAbs_(a, *this);
            sign_ = !sign_;
        }
    } else {
        addAbs_(*this, a);
    }
    normalize_();

    return *this;
}

BigInteger &BigInteger::operator*=(const BigInteger &a) {
    BigInteger res;
    res.numb_.resize(numb_.size() + a.numb_.size() + 2, 0);
    for (size_t i = 0; i < a.numb_.size(); ++i) {
        long long over = 0;
        for (size_t j = 0; j < numb_.size(); ++j) {
            long long tmp = static_cast<long long>(numb_[j]) * (a.numb_[i]) + over;
            over = tmp / static_cast<long long>(BASE_);
            tmp -= BASE_ * over;
            res.numb_[i + j] += tmp;
            if (res.numb_[i + j] >= BASE_) {
                over++;
                res.numb_[i + j] -= BASE_;
            }
        }
        int pos = i + numb_.size();
        while (over > 0) {
            res.numb_[pos] += over;
            over = 0;
            if (res.numb_[pos] >= BASE_) {
                over++;
                res.numb_[pos] -= BASE_;
            }
            pos++;
        }
    }
    res.deleteLeadingZeroes_();
    res.sign_ = !(sign_ ^ a.sign_);
    *this = res;
    normalize_();
    return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &b) {
    BigInteger res;
    BigInteger mod;
    div_(*this, b, res, mod);
    *this = res;
    normalize_();
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &b) {
    BigInteger res;
    BigInteger mod;
    div_(*this, b, res, mod);
    *this = mod;
    normalize_();
    return *this;
}

bool operator>(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) == 1;
}

bool operator<(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) == -1;
}

bool operator==(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) == 0;
}

bool operator!=(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) != 0;
}

bool operator<=(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) < 1;
}

bool operator>=(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) > -1;
}

BigInteger operator+(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c += b;
    return c;
}

BigInteger operator-(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c -= b;
    return c;
}

BigInteger operator*(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c *= b;
    return c;
}

BigInteger operator/(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c /= b;
    return c;
}

BigInteger operator%(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c %= b;
    return c;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &b) {
    return out << b.toString();
}

BigInteger &BigInteger::operator++() {
    return (*this) += 1;
}

BigInteger &BigInteger::operator--() {
    return (*this) -= 1;
}

BigInteger BigInteger::operator++(int) {
    BigInteger ans = (*this);
    ++(*this);
    return ans;
}

BigInteger BigInteger::operator-() const {
    BigInteger ans = *this;
    ans.sign_ ^= 1;
    return ans;
}

BigInteger BigInteger::operator--(int) {
    BigInteger ans = (*this);
    --(*this);
    return ans;
}

void BigInteger::normalize_() {
    if (compare(*this, 0, true) == 0) {
        sign_ = true;
    }
}

int BigInteger::getRes_(int a, int b) const {
    if (a < b) {
        return -1;
    }
    if (a > b) {
        return 1;
    }
    return 0;
}

std::istream &operator>>(std::istream &in, BigInteger &b) {
    std::string s;
    in >> s;

    return in;
}

BigInteger get_gcd(BigInteger a, BigInteger b) {
    while (b) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}

struct Rational {
public:
    Rational();

    Rational(long long a);

    Rational(int a);

    Rational(const BigInteger &a);

    std::string toString() const;

    explicit operator double();

    std::string asDecimal(size_t precision = 0) const;

    int compare(const Rational &a, bool byAbs) {
        return x_.compare(x_ * a.y_, a.x_ * y_, byAbs);
    }

    Rational &operator+=(const Rational &a);

    Rational &operator-=(const Rational &a);

    Rational &operator*=(const Rational &a);

    Rational &operator/=(const Rational &a);

    bool operator>(const Rational &a);

    bool operator<(const Rational &a);

    bool operator==(const Rational &a);

    bool operator!=(const Rational &a);

    bool operator<=(const Rational &a);

    bool operator>=(const Rational &a);

    friend Rational operator+(const Rational &a, const Rational &b);

    friend Rational operator-(const Rational &a, const Rational &b);

    friend Rational operator*(const Rational &a, const Rational &b);

    friend Rational operator/(const Rational &a, const Rational &b);

    Rational operator-() const;

private:
    BigInteger x_;
    BigInteger y_;

    void simplify_();
};

void Rational::simplify_() {
    BigInteger gcd = get_gcd(x_.abs(), y_.abs());
    x_ /= gcd;
    y_ /= gcd;

}

Rational::Rational() {
    x_ = 0;
    y_ = 1;
}

Rational::Rational(long long a) {
    x_ = a;
    y_ = 1;
}

Rational::Rational(int a) {
    x_ = a;
    y_ = 1;
}

Rational::              Rational(const BigInteger &a) {
    x_ = a;
    y_ = 1;
}

std::string Rational::toString() const {

    std::string ans;
    ans += x_.toString();
    if (y_ != 1 && x_ != 0) {
        ans += "/";
        ans += y_.toString();
    }
    return ans;
}


Rational::operator double() {
    std::string s = asDecimal(1000);
    double ans = std::stod(s);
    return ans;
}

Rational &Rational::operator+=(const Rational &a) {
    x_ = (x_ * a.y_) + (a.x_ * y_);
    y_ = y_ * a.y_;
    simplify_();
    return *this;
}

Rational &Rational::operator-=(const Rational &a) {
    x_ = (x_ * a.y_) - (a.x_ * y_);
    y_ = y_ * a.y_;
    simplify_();
    return *this;
}

Rational operator+(const Rational &a, const Rational &b) {
    Rational c = a;
    c += b;
    return c;
}

Rational operator-(const Rational &a, const Rational &b) {
    Rational c = a;
    c -= b;
    return c;
}

Rational operator*(const Rational &a, const Rational &b) {
    Rational c = a;
    c *= b;
    return c;
}

Rational operator/(const Rational &a, const Rational &b) {
    Rational c = a;
    c /= b;
    return c;
}

Rational Rational::operator-() const {
    Rational ans = *this;
    ans.x_ = -ans.x_;
    return ans;
}

Rational &Rational::operator/=(const Rational &a) {
    if (&a == this) return (*this) = 1;
    x_ *= a.y_;
    y_ *= a.x_;
    if (y_ < 0) {
        y_ = y_.abs();
        x_ = -x_;
    }
    simplify_();
    return *this;
}

Rational &Rational::operator*=(const Rational &a) {
    x_ *= a.x_;
    y_ *= a.y_;
    simplify_();
    return *this;
}

bool Rational::operator>(const Rational &a) {
    return (*this).compare(a, false) == 1;
}

bool Rational::operator<(const Rational &a) {
    return (*this).compare(a, false) == -1;
}

bool Rational::operator==(const Rational &a) {
    return (*this).compare(a, false) == 0;
}

bool Rational::operator!=(const Rational &a) {
    return (*this).compare(a, false) != 0;
}

bool Rational::operator<=(const Rational &a) {
    return (*this).compare(a, false) != 1;
}

bool Rational::operator>=(const Rational &a) {
    return (*this).compare(a, false) != -1;
}

std::string Rational::asDecimal(size_t precision) const {
    BigInteger mul = 1;
    for (size_t i = 0; i < precision; ++i) {
        mul *= 10;
    }
    BigInteger res = x_ * BigInteger(mul) / y_;
    std::string sign_string;
    if (res < 0) {
        res = -res;
        sign_string = "-";
    }
    std::string pre_ans = (res).toString();
    my_reverse(pre_ans.begin(), pre_ans.end());
    pre_ans.resize(std::max(pre_ans.size(), precision + 1), '0');
    my_reverse(pre_ans.begin(), pre_ans.end());
    std::string ans = sign_string + pre_ans.substr(0, pre_ans.size() - precision);
    if (precision != 0) {
        return ans + "." + pre_ans.substr(pre_ans.size() - precision, pre_ans.length());
    } else {
        return ans;
    }
}

#endif