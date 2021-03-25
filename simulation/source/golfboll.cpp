// Björn Sundin, 2021

#include <common.hpp>

//------------------------------------------

namespace constants {

// https://sv.wikipedia.org/wiki/Golfboll
constexpr auto ball_radius = 22.q_mm;
constexpr auto ball_mass = 45.q_g;

// https://en.wikipedia.org/wiki/Drag_coefficient
// 0.1 för mer turbulent, 0.47 för mer laminärt...
constexpr auto drag_coefficient = 0.2l;

constexpr auto cross_sectional_area = ball_radius*ball_radius*std::numbers::pi_v<long double>;

constexpr auto air_resistance = drag_coefficient*air_density*cross_sectional_area/2;

//------------------------------------------

constexpr auto start_speed = 40.q_m_per_s;
constexpr auto start_angle = degrees(45);

// Positivt är motsols.
constexpr auto revolutions_per_second = 4.l;
constexpr auto angular_velocity = revolutions_per_second*std::numbers::pi_v<long double>*2.l/1.q_s;
auto const magnus_force_coefficient = 2.l*std::numbers::pi_v<long double>*air_density*angular_velocity*units::pow<3>(ball_radius);

//------------------------------------------

// constexpr auto time_steps = std::array{1.q_ms, 5.q_ms, 10.q_ms, 50.q_ms};
constexpr auto time_steps = std::array{1.q_ms};

constexpr auto time = 6.q_s;

} // namespace constants

//------------------------------------------

struct Simulation 
{
    bool is_damped{};
    bool is_spinning{};
    bool is_euler_cromer{};
    
    Position position{};

    Velocity velocity{polar(constants::start_speed, constants::start_angle)};

    void step(units::physical::Time auto const time_step) 
    {
        auto const update_position = [=, this] { position += velocity*time_step; };
        
        if (!is_euler_cromer) {
            update_position();
        }

        if (is_spinning) {
            auto const magnus_force = constants::magnus_force_coefficient*Vec2{-velocity.y, velocity.x};
            velocity += magnus_force/constants::ball_mass*time_step;
        }
        if (is_damped) {
            auto const damping_force = constants::air_resistance*velocity.length_squared();
            auto const damping_acceleration = -velocity.normalized()*(damping_force/constants::ball_mass);
            velocity += damping_acceleration*time_step;
        }
        velocity += constants::g*time_step;

        if (is_euler_cromer) {
            update_position();
        }
    }
};

//------------------------------------------

decltype(auto) plot_simulation(Simulation simulation, units::physical::Time auto const time_step) 
{
    auto const iterations = calculate_iterations(time_step, constants::time);
    auto x = std::vector<long double>(iterations);
    auto y = std::vector<long double>(iterations);

    for (auto [a, b] : pairs(x, y)) 
    {
        a = simulation.position.x.count();
        b = simulation.position.y.count();
        simulation.step(time_step);
    }
    return matplot::plot(x, y)->line_width(3.f);
}

void plot_all() 
{
    matplot::colororder(std::vector<std::vector<float>>{
        {1.f, 0.4f, 0.4f}, //{0.4f, 0.f, 0.f},
        {0.4f, 1.f, 0.4f}, //{0.f, 0.4f, 0.f},
        {0.4f, 0.4f, 1.f}, //{0.f, 0.f, 0.4f}
    });

    for (auto const time_step : constants::time_steps) 
    {
        plot_simulation({.is_damped = false, .is_spinning = false, .is_euler_cromer = false}, time_step);
        // plot_simulation({.is_damped = false, .is_spinning = false, .is_euler_cromer = true}, time_step);
        plot_simulation({.is_damped = true, .is_spinning = false, .is_euler_cromer = false}, time_step);
        // plot_simulation({.is_damped = true, .is_spinning = false, .is_euler_cromer = true}, time_step);
        plot_simulation({.is_damped = true, .is_spinning = true, .is_euler_cromer = false}, time_step);
        // plot_simulation({.is_damped = true, .is_spinning = true, .is_euler_cromer = true}, time_step);
    }

    auto const legend = matplot::legend({"F_g", "F_g+F_D", "F_g+F_D+F_M"});
    legend->box(false);
    legend->location(matplot::legend::general_alignment::bottomright);
}

//------------------------------------------

int main() 
{
    matplot::gcf()->reactive_mode(false);
    matplot::gcf()->size(700, 400);
    matplot::hold(true);
    matplot::grid(true);
    
    matplot::xlabel("X-position [meter]");
    matplot::ylabel("Y-position [meter]");

    plot_all();

    matplot::axis(matplot::tight);

    matplot::save("results/golfboll.svg");
}
