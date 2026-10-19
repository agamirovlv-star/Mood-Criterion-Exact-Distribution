#include <bits/stdc++.h>
using namespace std;
using u128 = unsigned __int128;

// ======================================================================

// Вывод u128
ostream& operator<<(ostream& os, u128 v) {
    if (v == 0) { os << '0'; return os; }
    char buf[64]; int i = 0;
    while (v > 0) { buf[i++] = char('0' + (int)(v % 10)); v /= 10; }
    while (i--) os << buf[i];
    return os;
}

// =====================================================================

using Poly = vector<pair<long long, u128>>;

// dst += src, сдвинутый на a, с множителем mult (1 или 2)
static void add_shifted(Poly &dst, const Poly &src, long long a, int mult) {
    if (src.empty()) return;
    Poly res;
    res.reserve(dst.size() + src.size());
    size_t i = 0, k = 0;
    while (i < dst.size() || k < src.size()) {
        long long ds = (i < dst.size()) ? dst[i].first : LLONG_MAX;
        long long ss = (k < src.size()) ? src[k].first + a : LLONG_MAX;
        if (ds < ss) {
            res.push_back(dst[i++]);
        } else if (ss < ds) {
            u128 v = src[k].second;
            if (mult == 2) v <<= 1;
            res.emplace_back(ss, v);
            ++k;
        } else {
            u128 v = dst[i].second;
            u128 w = src[k].second;
            if (mult == 2) w <<= 1;
            res.emplace_back(ds, v + w);
            ++i; ++k;
        }
    }
    dst.swap(res);
}

// =====================================================================

int main() {
 
    ofstream fout("mood.out");
    if (!fout) return 1;

    auto start_time = chrono::high_resolution_clock::now();

    ifstream inp("mood.inp");
    if (!inp) { fout << "Cannot open mood.inp\n"; return 1; }

    string st;
    int k;
    inp >> st >> k >> st;
    vector<int> n_groups(k);
    for (int i = 0; i < k; ++i) inp >> n_groups[i];
    inp.close();

    int m=n_groups[0];
    int n=n_groups[1];

    int N = m + n;

    // a_i = (2i - N - 1)^2, i = 1..N
    vector<long long> a(N+1);
    for (int i = 1; i <= N; ++i) {
        long long t = 2LL*i - N - 1;
        a[i] = t*t;
    }

    // P[j] — разреженный полином
    vector<Poly> P(m+1);
    P[0].emplace_back(0, 1);

    // Обрабатываем пары (i, N+1-i), i = 1..N/2
    // Для нечётного N — центральный одиночный a[(N+1)/2] = 0
    int half = N / 2;
    for (int i = 1; i <= half; ++i) {
        long long ai = a[i];
        // пара: (1 + x^{ai} y)^2 = 1 + 2 x^{ai} y + x^{2ai} y^2
        int jmax = min(2*i, m);
        int jmin = max(0, m - (N - 2*i));  // не считаем то, что не успеет
        for (int j = jmax; j >= max(1, jmin); --j) {
            // P[j] += 2 * shift(P[j-1], ai)
            add_shifted(P[j], P[j-1], ai, 2);
            if (j >= 2) {
                // P[j] += shift(P[j-2], 2*ai)
                add_shifted(P[j], P[j-2], 2*ai, 1);
            }
        }
        // P[0] не меняется: (1 + ... )^2 даёт 1 при y^0
    }
    if (N % 2 == 1) {
        // центральный множитель (1 + x^0 y) = (1 + y), a = 0
        // P[j] += P[j-1] для всех j
        for (int j = min(half*2+1, m); j >= 1; --j) {
            add_shifted(P[j], P[j-1], 0, 1);
        }
    }


    const Poly &freq = P[m];
    u128 total = 0;
    for (auto &kv : freq) total += kv.second;

  
    fout << "Size ; " << freq.size() << "\n";
    fout << "m ; " << m << " ; n ; " << n << "\n";
    fout << "Freq ; Rang ; Prob ;\n";
    fout << setprecision(47);
       
     u128 acc = 0;
    for (auto &kv : freq) {
        acc += kv.second;
        long double prob = (long double)acc / (long double)total;
        long long rang = kv.first / 4;
        fout << rang << " ; " << (double)prob << "\n";
    }



    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    double time_elapsed = elapsed.count();
    fout << "Вычисление завершено за " << fixed << time_elapsed << " секунд\n";

   fout.close();
    return 0;
}