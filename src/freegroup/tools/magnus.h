#include <algorithm>
#include <cassert>
#include <cstdint>
#include <deque>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <optional>
#include <limits>

using Letter = int;
using Word = std::vector<Letter>;

template <class Vec> Vec Reciprocal(const Vec &word) {
    Vec result;
    result.reserve(word.size());
    for (std::size_t i = 0; i < word.size(); ++i) {
        result.push_back(-word[word.size() - i - 1]);
    }
    return result;
}

int cnt = 0;
int rw_cnt = 0;


template <class Vec>
bool IsTrivial(const Vec &word, const Vec &closure = {}) {
    // if (cnt++ % 5000 == 0) {
    //     std::cout << "Triv" << cnt << "\n";
    // }

    if (closure.empty()) {
        return word.empty();
    }

    std::size_t wlen = word.size();
    std::size_t clen = closure.size();
    if (closure.size() != 0 && word.size() > closure.size() && wlen > 2*clen) {
        return false;
    }

    auto inv = Reciprocal(closure);

    std::size_t maxOffset = 2*clen - wlen;

    for (std::size_t i = 0; i < wlen; ++i) {
        bool ok = true;
        for (std::size_t j = 0; j < wlen; ++j) {
            if (word[j] != closure[(i + j) % clen]) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return true;
        }

        ok = true;
        for (std::size_t j = 0; j < wlen; ++j) {
            if (word[j] != inv[(i + j) % clen]) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return true;
        }
    }

    return false;
}


template <class Vec, class El>
void ReduceModuloNormalClosureStep(Vec &reduced, El token,
                                   const Vec &closure = {}) {
    reduced.push_back(token);

    if (reduced.size() >= 2 && reduced[reduced.size() - 2] == -reduced.back()) {
        reduced.pop_back();
        reduced.pop_back();
    }

    if (!closure.empty() && reduced.size() >= closure.size()) {
        auto it_begin = reduced.end() - closure.size();
        Vec tail(it_begin, reduced.end());

        if (IsTrivial(tail, closure)) {
            reduced.erase(it_begin, reduced.end());
        }
    }
}

template <class Vec>
Vec ReduceModuloNormalClosure(const Vec &word, const Vec &closure = {}) {
    Vec reduced;

    for (auto token : word) {
        ReduceModuloNormalClosureStep(reduced, token, closure);
    }

    return reduced;
}

template <class Vec> Vec Normalize(const Vec &word) {
    return ReduceModuloNormalClosure(word);
}

inline std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return  x ^ (x >> 31);
}

inline void hash_combine(std::size_t& seed, std::size_t value) noexcept
{
    seed ^= splitmix64(value) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

class LetterWithSubscript {
public:
    LetterWithSubscript(Letter x, Word subs = {})
        : val_(x), subscript_(std::move(subs)), hash_(ComputeHash()) {}

    Letter Head() const { return val_; }

    std::size_t Hash() const noexcept { return hash_; }

    int Sign() const { return val_ > 0 ? 1 : -1; }

    LetterWithSubscript Abs() const { return Sign() == 1 ? *this : -*this; }

    const Word &Subscript() const { return subscript_; }

    Word Subscript() { return subscript_; }

    LetterWithSubscript operator-() const noexcept { return {static_cast<Letter>(-val_), subscript_}; }

    friend bool operator==(const LetterWithSubscript &a,
                            const LetterWithSubscript &b) noexcept {
        return a.val_ == b.val_ && a.subscript_ == b.subscript_;
    }

    friend bool operator!=(const LetterWithSubscript &a,
                            const LetterWithSubscript &b) noexcept {
        return !(a == b);
    }

    friend bool operator>(const LetterWithSubscript &a,
        const LetterWithSubscript &b) noexcept {
        if (a.val_ != b.val_) {
            return a.val_ > b.val_;
        } else {
            return a.subscript_ > b.subscript_;
        }
    }

    friend bool operator<(const LetterWithSubscript &a,
        const LetterWithSubscript &b) noexcept {
        if (a.val_ != b.val_) {
            return a.val_ < b.val_;
        } else {
            return a.subscript_ < b.subscript_;
        }
    }


    LetterWithSubscript operator+(Letter delta) const noexcept {
        auto new_sub = subscript_;
        if (new_sub.empty()) {
            new_sub.push_back(delta);
        } else {
            new_sub.back() += delta;
        }
        return { val_, std::move(new_sub) };
    }

    LetterWithSubscript operator-(Letter delta) const noexcept {
        auto new_sub = subscript_;
        if (new_sub.empty()) {
            new_sub.push_back(-delta);
        } else {
            new_sub.back() -= delta;
        }
        return { val_, std::move(new_sub) };
    }

    LetterWithSubscript RemoveSubscript() const {
        auto new_sub = subscript_;
        if (!new_sub.empty()) new_sub.pop_back();
        return LetterWithSubscript{ val_, new_sub };
    }

    LetterWithSubscript AddSubscript(Letter sub) const {
        auto new_sub = subscript_;
        new_sub.push_back(sub);
        return LetterWithSubscript{ val_, new_sub };
    }

private:
    std::size_t ComputeHash() const noexcept {
        std::uint64_t h = splitmix64(static_cast<std::uint64_t>(val_));   // base letter **with sign**

        h ^= splitmix64(static_cast<std::uint64_t>(subscript_.size()));   // vector length

        for (Letter s : subscript_) {
            h ^= splitmix64(static_cast<std::uint64_t>(s));               // every component
        }
        return static_cast<std::size_t>(h);
    }

    Letter val_;
    Word subscript_;
    std::size_t hash_;
};

using SubWord = std::vector<LetterWithSubscript>;

std::size_t CountOccurrences(Letter x, const Word &word) {
    std::size_t count = 0;
    for (Letter l : word) {
        if (std::abs(l) == x) {
            count++;
        }
    }
    return count;
}

std::size_t CountOccurrences(const LetterWithSubscript& x, const SubWord& word) {
    std::size_t count = 0;
    for (LetterWithSubscript l : word) {
        if (l.Abs() == x) {
            count++;
        }
    }
    return count;
}

int CountExponent(Letter x, const Word& word) {
    int count = 0;
    for (Letter l : word) {
        if (l == x) {
            count++;
        } else if (-l == x) {
            count--;
        }
    }
    return count;
}

int CountExponent(const LetterWithSubscript& x, const SubWord& word) {
    int count = 0;
    for (LetterWithSubscript l : word) {
        if (l == x) {
            count++;
        } else if (l == -x) {
            count--;
        }
    }
    return count;
}

int Sign(long long x) { return x > 0 ? 1 : -1; }

namespace std {
    template <> struct hash<LetterWithSubscript>
    {
        std::size_t operator()(const LetterWithSubscript& l) const noexcept
        {
            return l.Hash();
        }
    };

    // ----- SubWord --------------------------------------------------------------
    template <> struct hash<SubWord>
    {
        std::size_t operator()(const SubWord& w) const noexcept
        {
            std::size_t seed = w.size();         // distinguish [a,b] from [b,a]
            for (const auto& letter : w)
                hash_combine(seed, std::hash<LetterWithSubscript>{}(letter));
            return seed;
        }
    };
} // namespace std

std::pair<std::vector<SubWord>, std::vector<SubWord>>
SplitByT(const SubWord &word,
         const std::set<LetterWithSubscript> &setT) {
    if (word.empty()) {
        return {{{}}, {}};
    }

    std::vector<SubWord> ts({{}});
    std::vector<SubWord> splits({});

    bool prev_t = true;

    for (const LetterWithSubscript &v : word) {
        bool in_t = setT.count(v.Abs()) > 0;
        if (in_t && prev_t) {
            ts.back().push_back(v);
        } else if (in_t) {
            ts.push_back({v});
        } else if (prev_t) {
            splits.push_back({v});
        } else {
            splits.back().push_back(v);
        }
        prev_t = in_t;
    }

    if (setT.count(word.back().Abs()) == 0) {
        ts.push_back({});
    }

    return {ts, splits};
}

SubWord CyclicallyReduce(const SubWord &word) {
    std::size_t n = word.size();
    std::deque<std::pair<std::size_t, LetterWithSubscript>> stack;
    SubWord opt;

    for (std::size_t i = 0; i < 2 * n; ++i) {
        const LetterWithSubscript &x = word[i % word.size()];
        if (stack.empty()) {
            stack.push_back({i, x});
        } else {
            std::size_t j = stack.back().first;
            LetterWithSubscript &y = stack.back().second;

            if (y == -x) {
                stack.pop_back();
            } else {
                stack.push_back({i, x});
            }
        }

        if (stack.empty() && i >= n) {
            return {};
        }

        if (!stack.empty()) {
            std::size_t j = stack.front().first;
            LetterWithSubscript &y = stack.back().second;

            if (i >= j && i - j >= n) {
                stack.pop_front();
            }

            if (i >= n) {
                if (opt.empty() || stack.size() <= opt.size()) {
                    opt.clear();
                    for (const auto &elem : stack) {
                        opt.push_back(elem.second);
                    }
                }
            }
        }
    }

    return opt;
}

template <class El>
void AddNTimes(const El &T, std::ptrdiff_t n, std::vector<El> &source) {
    if (n == 0) {
        return;
    }

    auto count = static_cast<std::size_t>(std::llabs(n));
    source.reserve(source.size() + count);

    if (n > 0) {
        source.insert(source.end(), count, T);
    } else {
        El neg = -T;
        source.insert(source.end(), count, neg);
    }
}

template <class El> SubWord RemoveSubscript(const SubWord &word, const El &T) {
    SubWord w;
    auto neg = -T;

    for (const auto &v : word) {
        AddNTimes(T, v.Subscript().back(), w);
        w.push_back(v.RemoveSubscript());
        AddNTimes(neg, v.Subscript().back(), w);
    }

    return w;
}

SubWord Psi(const SubWord &word, const LetterWithSubscript &t,
            const LetterWithSubscript &x, int alpha, int beta) {
    SubWord w;

    for (const auto &v : word) {
        if (v.Abs() == x.Abs()) {
            if (v.Sign() > 0) {
                w.push_back(x);
                AddNTimes(t, -alpha, w);
            } else {
                AddNTimes(t, alpha, w);
                w.push_back(-x);
            }
        } else if (v.Abs() == t.Abs()) {
            AddNTimes(v, beta, w);
        } else {
            w.push_back(v);
        }
    }

    return w;
}

SubWord PsiPseudoInverse(const SubWord &word, const LetterWithSubscript &t,
                         const LetterWithSubscript &x, int alpha, int beta) {
    SubWord w;

    for (const auto &v : word) {
        if (!(v.Abs() == x.Abs())) {
            w.push_back(v);
        } else {
            if (v.Sign() > 0) {
                w.push_back(x);
                AddNTimes(t, alpha, w);
            } else {
                AddNTimes(t, -alpha, w);
                w.push_back(-x);
            }
        }
    }

    return w;
}

template <typename T> std::size_t CountUnique(const std::vector<T> &v) {
    std::set<T> seen;
    // seen.reserve(v.size());
    for (auto const &x : v) {
        seen.insert(x);
    }
    return seen.size();
}

std::pair<int, SubWord>
PropagateT(int cumsum, SubWord &word, const SubWord &relator,
           const LetterWithSubscript &x, Letter a, Letter b,
           std::map<SubWord, SubWord>& cache_a,
           std::map<SubWord, SubWord>& cache_b);

SubWord ReduceWordProblem(const SubWord &w, const SubWord &rel,
                          std::set<LetterWithSubscript>& T) {
    // if (rw_cnt++ % 100 == 0) {
    //     std::cout << "RW:" << rw_cnt << "\n";
    // }
    SubWord relator = CyclicallyReduce(rel);
    SubWord word = Normalize(w);
    std::map<LetterWithSubscript, std::vector<std::size_t>> letter_to_pos;

    for (std::size_t i = 0; i < relator.size(); ++i) {
        letter_to_pos[relator[i].Abs()].push_back(i);
    }

    if (word.empty() || std::all_of(word.begin(), word.end(),
                   [&](auto &v){ return T.count(v.Abs()) > 0; })) {
        return word;
    }

    bool flag = true;

    for (auto v : letter_to_pos) {
        if (T.count(v.first) == 0) {
            flag = false;
            break;
        }
    }

    if (flag) {
        auto ret = SplitByT(word, T);
        auto &flags = ret.first;
        auto &splits = ret.second;

        SubWord reduced_by_t;
        reduced_by_t.reserve(word.size());

        for (std::size_t i = 0; i < std::min(flags.size(), splits.size()); ++i) {
            if (!flags[i].empty()) {
                std::set<LetterWithSubscript> emptyT;
                auto mid_res = ReduceWordProblem(splits[i], relator, emptyT);
                reduced_by_t.insert(reduced_by_t.end(), mid_res.begin(), mid_res.end());
            } else {
                reduced_by_t.insert(reduced_by_t.end(), splits[i].begin(), splits[i].end());
            }
        }

        return Normalize(std::move(reduced_by_t));
    }

    if (CountUnique(relator) == 1) {
        return ReduceModuloNormalClosure(word, relator);
    }

    SubWord appropriate_ts;

    for (const LetterWithSubscript &v : relator) {
        if (CountExponent(v, relator) == 0 && CountOccurrences(v, relator) > 0) {
            appropriate_ts.push_back(v.Abs());
        }
    }

    if (!appropriate_ts.empty()) {
        LetterWithSubscript t = appropriate_ts[0];
        std::size_t x_idx;

        if (T.count(t) > 0) {
            for (std::size_t i = 0; i < relator.size(); ++i) {
                if (T.count(relator[i].Abs()) == 0) {
                    x_idx = i;
                    break;
                }
            }
        } else {
            for (std::size_t i = 0; i < relator.size(); ++i) {
                if (relator[i].Abs() != t) {
                    x_idx = i;
                    break;
                }
            }
        }

        std::rotate(relator.begin(), relator.begin() + x_idx, relator.end());
        LetterWithSubscript x = relator[0].Abs();

        assert(x != t);

        SubWord r_prime;
        int carry = 0;

        auto ret = SplitByT(relator, {t});
        auto &t_splits = ret.first;
        auto &splits = ret.second;

        for (std::size_t i = 0; i < std::min(t_splits.size(), splits.size()); ++i) {
            if (!t_splits[i].empty() && t_splits[i][0].Sign() > 0) {
                carry += static_cast<int>(t_splits[i].size());
            } else {
                carry -= static_cast<int>(t_splits[i].size());
            }

            SubWord cur_vec;
            cur_vec.reserve(splits[i].size());

            for (const LetterWithSubscript &elem : splits[i]) {
                cur_vec.push_back(elem.AddSubscript(static_cast<Letter>(carry)));
            }

            r_prime.insert(r_prime.end(), cur_vec.begin(), cur_vec.end());
        }

        assert(Normalize(RemoveSubscript(r_prime, t)) == relator);
        Letter a = std::numeric_limits<Letter>::max();
        Letter b = std::numeric_limits<Letter>::min();

        for (const LetterWithSubscript &v : r_prime) {
            if (v.Abs().RemoveSubscript() == x) {
                a = std::min(a, v.Subscript().back());
                b = std::max(b, v.Subscript().back());
            }
        }

        assert(a <= 0 && b >= 0);

        ret = SplitByT(word, {t});
        t_splits = ret.first;
        splits = ret.second;

        std::vector<int> cumsums;
        cumsums.reserve(t_splits.size());

        for (const auto &w : t_splits) {
            if (!w.empty() && w[0].Sign() > 0) {
                cumsums.push_back(static_cast<int>(w.size()));
            } else {
                cumsums.push_back(-static_cast<int>(w.size()));
            }
        }

        for (std::size_t i = 0; i < splits.size(); ++i) {
            for (std::size_t j = 0; j < splits[i].size(); ++j) {
                splits[i][j] = splits[i][j].AddSubscript(0);
            }
        }

        std::size_t iteration = 0;

        std::map<SubWord, SubWord> cache_a;
        std::map<SubWord, SubWord> cache_b;

        while (true) {
            int carry = 0;
            bool flag = false;

            std::vector<int> cumsums_;
            std::vector<SubWord> splits_;

            for (std::size_t i = 0; i < std::min(cumsums.size(), splits.size()); ++i) {
                int cumsum = cumsums[i];
                auto ret2 = PropagateT(carry + cumsum, splits[i], r_prime, x, a, b, cache_a, cache_b);
                int carry_ = ret2.first;
                SubWord w = std::move(ret2.second);

                if (!cumsums_.empty() && carry + cumsum == carry_) {
                    splits_.back().insert(splits_.back().end(), w.begin(), w.end());
                } else {
                    cumsums_.push_back(carry + cumsum - carry_);
                    splits_.push_back(w);
                }

                carry = carry_;
                flag = flag || (carry_ != 0);
            }

            cumsums_.push_back(cumsums.size() > splits.size() ? carry + cumsums.back() : 0);
            cumsums = cumsums_;
            splits = splits_;

            if (!flag) {
                break;
            }

            iteration++;
        }

        std::set<LetterWithSubscript> T_prime;
        if (T.count(t.Abs()) > 0) {
            for (const SubWord& spl : splits) {
                for (const LetterWithSubscript& v : spl) {
                    if (T.count(v.RemoveSubscript().Abs()) > 0) {
                        T_prime.insert(v.Abs());
                    }
                }
            }

            for (const LetterWithSubscript& v : r_prime) {
                if (T.count(v.RemoveSubscript().Abs()) > 0) {
                    T_prime.insert(v.Abs());
                }
            }
        } else {
            for (const LetterWithSubscript& v : T) {
                T_prime.insert(v.AddSubscript(0));
            }
        }

        for (std::size_t i = 0; i < splits.size(); ++i) {
            splits[i] = ReduceWordProblem(splits[i], r_prime, T_prime);
        }

        SubWord result;

        for (std::size_t i = 0; i < std::min(cumsums.size(), splits.size()); ++i) {
            AddNTimes(t, cumsums[i], result);
            auto tmp = RemoveSubscript(splits[i], t);
            result.insert(result.end(), tmp.begin(), tmp.end());
        }

        AddNTimes(t, cumsums.back(), result);
        result = Normalize(result);

        return result;
    }

    LetterWithSubscript t(0);
    LetterWithSubscript x(0);

    std::vector<LetterWithSubscript> not_in_T;
    not_in_T.reserve(letter_to_pos.size());
    for (const LetterWithSubscript& le : relator) {
        if (T.find(le.Abs()) == T.end()) {
            not_in_T.push_back(le.Abs());
        }
    }

    if (not_in_T.size() >= 2) {
        t = not_in_T[0];
        x = not_in_T[1];
    } else {
        for (const LetterWithSubscript& le : relator) {
            if (T.count(le.Abs())) {
                t = le.Abs();
                break;
            }
        }

        for (const LetterWithSubscript& le : relator) {
            if (!(le.Abs() == t) && T.count(le.Abs()) == 0) {
                x = le.Abs();
                break;
            }
        }
    }

    int alpha = CountExponent(t, relator);
    int beta = CountExponent(x, relator);

    SubWord r_prime = Psi(relator, t, x, alpha, beta);
    SubWord w_prime = Psi(word, t, x, alpha, beta);

    SubWord result = ReduceWordProblem(w_prime, r_prime, T);

    if (T.count(t) == 0) {
        result = Normalize(result);

        return result;
    }

    result = Normalize(PsiPseudoInverse(result, t, x, alpha, beta));

    SubWord ws;
    auto ret = SplitByT(result, {t});
    auto &ts = ret.first;
    auto &splits = ret.second;
    std::size_t total = ts.size() + splits.size();

    for (std::size_t i = 0; i < total; ++i) {
        // Processing splits
        std::size_t real_index = i / 2;
        if (i % 2 != 0) {
            if (real_index >= splits.size()) {
                continue;
            }
            ws.insert(ws.end(), splits[real_index].begin(), splits[real_index].end());
        } else {
            if (real_index >= ts.size()) {
                continue;
            }
            if (ts[real_index].empty()) {
                continue;
            }

            if (ts[real_index].size() % std::abs(beta) == 0) {
                AddNTimes(ts[real_index][0], Sign(beta) * ts[real_index].size() / std::abs(beta), ws);
            } else {
                return {};
            }
        }
    }

    return ws;
}

bool ContainsAbs(const SubWord& vec, LetterWithSubscript val) {
    for (const LetterWithSubscript& l : vec) {
        if (l.Abs() == val) {
            return true;
        }
    }
    return false;
}

std::pair<int, SubWord>
PropagateT(int cumsum, SubWord& word, const SubWord &relator,
           const LetterWithSubscript &x, Letter a, Letter b,
           std::map<SubWord, SubWord>& cache_a,
           std::map<SubWord, SubWord>& cache_b) {
    LetterWithSubscript x_a = x.AddSubscript(a);
    LetterWithSubscript x_b = x.AddSubscript(b);

    std::set<LetterWithSubscript> setA;
    std::set<LetterWithSubscript> setB;

    for (std::size_t i = 0; i < word.size(); ++i) {
        auto ab = word[i].Abs();
        if (ab != x_a) {
            setA.insert(ab);
        }
        if (ab != x_b) {
            setB.insert(ab);
        }
    }

    for (std::size_t i = 0; i < relator.size(); ++i) {
        auto ab = relator[i].Abs();
        if (ab != x_a) {
            setA.insert(ab);
        }
        if (ab != x_b) {
            setB.insert(ab);
        }
    }

    int carry = 0;
    int iteration = 0;

    while (cumsum != 0) {
        for (std::size_t i = 0; i < word.size(); ++i) {
            auto ab = word[i].Abs();
            if (ab != x_a) {
                setA.insert(ab);
            }
            if (ab != x_b) {
                setB.insert(ab);
            }
        }

        if (cumsum > 0 && ContainsAbs(word, x_b)) {
            auto it = cache_b.find(word);
            if (it == cache_b.end()) {
                it = cache_b.emplace(word,
                                     ReduceWordProblem(word, relator, setB)).first;
            }
            SubWord word_prime = it->second;

            if (!word_prime.empty()) {
                bool in_b = true;
                for (const LetterWithSubscript& le : word_prime) {
                    if (setB.find(le.Abs()) == setB.end()) {
                        in_b = false;
                        break;
                    }
                }

                if (in_b) {
                    word = word_prime;
                } else {
                    break;
                }
            } else {
                break;
            }
        } else if (cumsum > 0) {
            for (std::size_t i = 0; i < word.size(); ++i) {
                word[i] = word[i] + 1;
            }
            cumsum--;
            carry++;
        } else if (cumsum < 0 && ContainsAbs(word, x_a)) {
            auto it = cache_a.find(word);
            if (it == cache_a.end()) {
                it = cache_a.emplace(word,
                                     ReduceWordProblem(word, relator, setA)).first;
            }
            SubWord word_prime = it->second;

            if (!word_prime.empty()) {
                bool in_a = true;
                for (const LetterWithSubscript& le : word_prime) {
                    if (setA.find(le.Abs()) == setA.end()) {
                        in_a = false;
                        break;
                    }
                }

                if (in_a) {
                    word = word_prime;
                } else {
                    break;
                }
            } else {
                break;
            }
        } else if (cumsum < 0) {
            for (std::size_t i = 0; i < word.size(); ++i) {
                word[i] = word[i] - 1;
            }
            cumsum++;
            carry--;
        }

        iteration++;
    }

    return { carry, word };
}

SubWord MagnusReduceModuloNormalClosure(const Word& w, const Word& r, std::set<Letter> t = {}) {
    SubWord word;
    SubWord relator;
    std::set<LetterWithSubscript> T;

    word.reserve(w.size());
    relator.reserve(r.size());

    for (std::size_t i = 0; i < w.size(); ++i) {
        word.push_back(LetterWithSubscript(w[i]));
    }

    for (std::size_t i = 0; i < r.size(); ++i) {
        relator.push_back(LetterWithSubscript(r[i]));
    }

    for (const Letter& l : t) {
        T.insert(LetterWithSubscript(l).Abs());
    }

    return ReduceWordProblem(word, relator, T);
}

bool MagnusIsFromNormalClosure(const Word& w, const Word& r, std::set<Letter> t = {}) {
    std::set<LetterWithSubscript> T;
    for (const Letter& l : t) {
        T.insert(LetterWithSubscript(l).Abs());
    }

    SubWord word = MagnusReduceModuloNormalClosure(w, r, t);

    if (word.empty()) {
        return true;
    }

    bool flag = true;

    for (const LetterWithSubscript& le : word) {
        if (T.count(le.Abs()) == 0) {
            return false;
        }
    }
    return true;
}
