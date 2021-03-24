#pragma once

#include <numbers>
#include <chrono>
#include <thread>

#include <matplot/matplot.h>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <units/physical/si.h>
#include <units/format.h>
#include <units/math.h>

//------------------------------------------

using namespace std::chrono_literals;

namespace si = units::physical::si;
using namespace si::literals;

//------------------------------------------

constexpr auto degrees(long double degrees) {
    return si::angle<si::radian, long double>{degrees*std::numbers::pi_v<long double>/180.l};
}

constexpr int calculate_iterations(units::physical::Time auto const time_step, units::physical::Time auto const time = 1.q_s) {
    return units::quantity_cast<units::dimensionless<units::unitless, int>>(time/time_step).count();
}

//------------------------------------------

constexpr auto range(std::integral auto const n) {
    return std::views::iota(0, n);
}
template<std::integral T>
constexpr auto range(T const start, T const end) {
    return std::views::iota(start, end);
}

template<std::ranges::random_access_range A, std::ranges::random_access_range B>
auto pairs(A& a, B& b) {
    return std::views::transform(a, [b_iterator = std::ranges::begin(b)](auto& a_element) mutable {
        return std::tie(a_element, *(b_iterator++));
    });
}

//------------------------------------------

template<units::Quantity T>
struct Vec2 {
    T x, y;

    constexpr T length_squared() const {
        return x*x + y*y;
    }
    T length() const {
        return T{std::hypot(x.count(), y.count())};
    }
    auto normalized() const {
        auto const len = length();
        return Vec2<decltype(x/len)>{x/len, y/len};
    }

    Vec2 rotated(si::angle<si::radian> const angle) {
        return Vec2{
            std::cos(angle.count())*x - std::sin(angle.count())*y,
            std::sin(anlge.count())*x + std::cos(angle.count())*y
        };
    }

    constexpr Vec2 operator-() const {
        return Vec2{-x, -y};
    }

    template<typename U>
    constexpr Vec2& operator+=(Vec2<U> const other) {
        x += units::quantity_cast<typename T::unit>(other.x);
        y += units::quantity_cast<typename T::unit>(other.y);
        return *this;
    }
};
template<typename T>
Vec2(T, T) -> Vec2<T>;

template<units::Quantity A, units::Quantity B>
constexpr auto operator+(Vec2<A> const lhs, Vec2<B> const rhs) {
    return Vec2{lhs.x + rhs.x, lhs.y + rhs.y};
}

template<units::Quantity A, units::Quantity B>
constexpr auto operator*(A const lhs, Vec2<B> const rhs) {
    return Vec2{lhs*rhs.x, lhs*rhs.y};
}
template<units::Quantity A, units::Quantity B>
constexpr auto operator*(Vec2<A> const lhs, B const rhs) {
    return Vec2{lhs.x*rhs, lhs.y*rhs};
}
template<units::Quantity A, units::Quantity B>
constexpr auto operator/(Vec2<A> const lhs, B const rhs) {
    return Vec2{lhs.x/rhs, lhs.y/rhs};
}

auto polar(units::Quantity auto magnitude, units::physical::Angle auto angle) {
    return Vec2{magnitude*std::cos(angle.count()), magnitude*std::sin(angle.count())};
}

//------------------------------------------

using PositionScalar = si::length<si::metre, long double>;
using VelocityScalar = si::speed<si::metre_per_second, long double>; 
using AccelerationScalar = si::acceleration<si::metre_per_second_sq, long double>;

using Position = Vec2<PositionScalar>;
using Velocity = Vec2<VelocityScalar>;
using Acceleration = Vec2<AccelerationScalar>;

//------------------------------------------

namespace constants {

// https://en.wikipedia.org/wiki/Density_of_air
constexpr auto air_density = 1.225q_kg_per_m3;

constexpr auto g = Vec2{0.q_m_per_s2, -units::physical::si::si2019::standard_gravity<long double>};

} // namespace constants

