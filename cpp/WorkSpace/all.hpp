#ifndef _ALL_H
#define _ALL_H

// C
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cuchar>
#include <cwchar>
#include <cwctype>

// C++ I/O
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <streambuf>

// 容器
#include <array>
#include <bitset>
#include <deque>
#include <forward_list>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// 算法与实用工具
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <random>
#include <ratio>
#include <tuple>
#include <utility>
#include <variant>
#include <optional>
#include <any>
#include <bit>
#include <compare>
#include <execution>
#include <ranges>
#include <span>
#include <string>
#include <string_view>

// 多线程与原子操作
#include <atomic>
#include <condition_variable>
#include <future>
#include <mutex>
#include <shared_mutex>
#include <thread>

// 时间与时区
#include <chrono>

// 类型信息与异常
#include <exception>
#include <stdexcept>
#include <type_traits>
#include <typeindex>
#include <typeinfo>

// 文件系统
#include <filesystem>

// 内存管理
#include <memory>
#include <memory_resource>
#include <new>
#include <scoped_allocator>

// 数值库
#include <complex>
#include <valarray>
#include <limits>

// 正则
#include <regex>

// 随机与数学
#include <random>

// 输入输出
#include <ios>
#include <iosfwd>
#include <ostream>
#include <istream>

// C++20/23 扩展（可选）
#ifdef __cpp_lib_format
#include <format>
#endif
#ifdef __cpp_lib_coroutine
#include <coroutine>
#endif
#ifdef __cpp_lib_semaphore
#include <semaphore>
#endif
#ifdef __cpp_lib_barrier
#include <barrier>
#endif
#ifdef __cpp_lib_latch
#include <latch>
#endif

using namespace std;

#endif // _ALL_H
