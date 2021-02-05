//
// Created by vahna on 28.10.2019.
//

#ifndef BIGINTEGER_PERMUTATION_H
#define BIGINTEGER_PERMUTATION_H

template<class curr>
void swap(curr &a, curr &b) {
    curr tmp = a;
    a = b;
    b = tmp;
}

template<class curr>
void reverse(curr first, curr last) {
    while ((first != last) && (first != --last)) {
        swap(*first++, *last);
    }
}

int min(const int a, const int b) {
    return a < b ? a : b;
}

class Permutation {
public:
    explicit Permutation(unsigned int n);

    Permutation(unsigned int n, const int *from);

    Permutation(const Permutation &from);

    Permutation &operator=(const Permutation &x);

    Permutation inverse() const;

    void copyArr(const int *from, unsigned int new_length);

    Permutation &operator++();

    Permutation &operator--();

    Permutation operator++(int);

    Permutation operator--(int);
    void copyArrSafe(int *to,const int *from) const;

    friend Permutation operator*(const Permutation &a, const Permutation &b);

    Permutation &operator*=(const Permutation &a);

    int compare(const Permutation &a) const {
        for (unsigned int i = 0; i < length_; i++) {
            if (all_[i] != a.all_[i]) {
                return all_[i] < a.all_[i] ? -1 : 1;
            }
        }
        return 0;
    }

    bool operator<(const Permutation &a);

    bool operator>(const Permutation &a);

    bool operator==(const Permutation &a);

    bool operator!=(const Permutation &a);

    bool operator>=(const Permutation &a);

    bool operator<=(const Permutation &a);

    int operator[](int i);

    void operator()(int *to) const;

    Permutation next();

    Permutation previous();

    ~Permutation();

private:
    int *all_;
    unsigned int length_;

    void multiplyMinus_();

};

Permutation::~Permutation() {
    delete[] all_;
}

Permutation Permutation::next() {
    Permutation ans(*this);
    ++ans;
    return Permutation(ans);
}

Permutation Permutation::previous() {
    Permutation ans(*this);
    --ans;
    return Permutation(ans);
}

int Permutation::operator[](int i) {
    return all_[i];
}

void Permutation::operator()(int *to) const {
    int *ans = new int[length_];
    for (unsigned int i = 0; i < length_; ++i) {
        ans[all_[i]] = to[i];
    }

    copyArrSafe(to,ans);
    delete[]ans;
}

Permutation Permutation::operator--(int) {
    Permutation ans(*this);
    --(*this);
    return ans;
}

Permutation Permutation::operator++(int) {
    const Permutation ans(*this);
    ++(*this);
    return ans;
}

Permutation &Permutation::operator--() {
    multiplyMinus_();
    ++*(this);
    multiplyMinus_();
    return *this;
}

Permutation &Permutation::operator++() {
    int min = 0;
    for (int i = static_cast<int>(length_) - 2; i >= 0; --i) {

        if (all_[i] < all_[i + 1]) {
            min = i + 1;
            for (unsigned int j = i + 1; j < length_; ++j) {
                if (all_[j] < all_[min] && all_[j] > all_[i]) {
                    min = static_cast<int>(j);
                }
            }
            swap(all_[i], all_[min]);
            reverse(all_ + i + 1, all_ + length_);
            return *this;
        }

    }
    return *this;

}

Permutation &Permutation::operator*=(const Permutation &a) {
    Permutation ans(length_);
    for (unsigned int i = 0; i < a.length_; ++i) {
        ans.all_[i] = all_[a.all_[i]];
    }
    *(this) = ans;
    return *(this);
}

Permutation Permutation::inverse() const {
    Permutation ans(length_);
    for (unsigned int i = 0; i < length_; ++i) {
        ans.all_[all_[i]] = static_cast<int>(i);
    }
    return ans;
}

Permutation::Permutation(unsigned int n) {
    length_ = n;
    all_ = new int[length_];
    for (unsigned int i = 0; i < length_; ++i) {
        all_[i] = static_cast<int>(i);
    }
}

Permutation::Permutation(unsigned int n, const int *from) {
    copyArr(from, n);
}

Permutation::Permutation(const Permutation &from) {
    if (this == &from) {
        length_ = 0;
        all_ = nullptr;
        return;
    }
    copyArr(from.all_, from.length_);
}

Permutation &Permutation::operator=(const Permutation &x) {
    if (this == &x) {
        return *this;
    }
    length_ = x.length_;
    delete[]all_;
    copyArr(x.all_,x.length_);
    return *this;

}

Permutation operator*(const Permutation &a, const Permutation &b) {
    Permutation ans(a);
    ans *= b;
    return ans;
}

void Permutation::multiplyMinus_() {
    for (unsigned int i = 0; i < length_; ++i) {
        all_[i] = -all_[i];
    }
}

void Permutation::copyArr(const int *from, unsigned int new_length) {
    length_ = new_length;
    all_ = new int[length_];
    copyArrSafe(all_,from);
}

bool Permutation::operator<(const Permutation &a) {
    return (*this).compare(a) == -1;
}

bool Permutation::operator>(const Permutation &a) {
    return (*this).compare(a) == 1;
}


bool Permutation::operator==(const Permutation &a) {
    return (*this).compare(a) == 0;
}

bool Permutation::operator!=(const Permutation &a) {
    return (*this).compare(a) != 1;
}

bool Permutation::operator>=(const Permutation &a) {
    return (*this).compare(a) != -1;
}

bool Permutation::operator<=(const Permutation &a) {
    return (*this).compare(a) != 1;
}

void Permutation::copyArrSafe( int *to,const int *from) const{
    for(unsigned int i = 0;i < length_;++i){
        to[i] = from[i];
    }
}

#endif //BIGINTEGER_PERMUTATION_H
