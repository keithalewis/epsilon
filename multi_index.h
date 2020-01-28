// multi_index.h - multiple index values
#pragma once
#include <algorithm>
#include <vector>
#include <map>

namespace fms {

	using multi_index = std::vector<int>;
	using multinomial = std::map<multi_index, double>;
}

inline bool operator==(const fms::multi_index& a, const fms::multi_index& b)
{
	return std::equal(a.begin(), a.end(), b.begin(), b.end());
}
inline bool operator!=(const fms::multi_index& a, const fms::multi_index& b)
{
	return !operator==(a, b);
}

inline fms::multi_index operator*(const fms::multi_index& a, const fms::multi_index& b)
{
	fms::multi_index c;

	if (a.size() > b.size()) {
		c = a;
		for (size_t i = 0; i < b.size(); ++i) {
			c[i] += b[i];
		}
	}
	else {
		c = b;
		for (size_t i = 0; i < a.size(); ++i) {
			c[i] += a[i];
		}
	}

	return c;
}
inline fms::multi_index operator/(const fms::multi_index& a, const fms::multi_index& b)
{
	size_t m = std::min(a.size(), b.size());
	size_t n = std::max(a.size(), b.size());
	fms::multi_index c(n);

	for (size_t i = 0; i < m; ++i) {
		c[i] = a[i] - b[i];
	}
	if (m == a.size()) {
		for (size_t i = m; i < n; ++i) {
			c[i] = -b[i];
		}
	}
	else {
		for (size_t i = m; i < n; ++i) {
			c[i] = a[i];
		}
	}

	return c;
}

inline fms::multinomial operator+(const fms::multinomial& a, const fms::multinomial& b)
{
	fms::multinomial c;

	for (const auto& [ka, va] : a) {
		for (const auto& [kb, vb] : b) {
			if (ka == kb) {
				c[ka] = va + vb;
			}
		}
	}

	return c;
}
inline fms::multinomial operator-(const fms::multinomial& a, const fms::multinomial& b)
{
	fms::multinomial c;

	for (const auto& [ka, va] : a) {
		for (const auto& [kb, vb] : b) {
			if (ka == kb) {
				c[ka] = va - vb;
			}
		}
	}

	return c;
}
inline fms::multinomial operator*(const fms::multinomial& a, const fms::multinomial& b)
{
	fms::multinomial c;

	for (const auto& [ka, va] : a) {
		for (const auto& [kb, vb] : b) {
			c[ka*kb] += va * vb;
		}
	}

	return c;
}