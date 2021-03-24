#include <common.hpp>

//------------------------------------------

namespace constants {
    
// https://sv.wikipedia.org/wiki/Bordtennis
constexpr auto ball_radius = 4q_cm;
constexpr auto ball_mass = 2.7q_g;

// https://en.wikipedia.org/wiki/Drag_coefficient
// 0.1 för mer turbulent, 0.47 för mer laminärt...
constexpr auto drag_coefficient = 0.3l;

constexpr auto cross_sectional_area = ball_radius*ball_radius*std::numbers::pi_v<long double>;

constexpr auto air_resistance = drag_coefficient*air_density*cross_sectional_area/2;

//------------------------------------------

constexpr auto start_speed = 40q_m_per_s;
constexpr auto start_angle = degrees(35);

//------------------------------------------

// constexpr auto time_steps = std::array{1.q_ms, 5.q_ms, 10.q_ms, 50.q_ms};
constexpr auto time_steps = std::array{0.1q_ms, 1.q_ms, 5.q_ms, 10.q_ms};

constexpr auto time = 2.q_s;

} // namespace constants

//------------------------------------------

struct Simulation 
{
    bool is_damped{};
    bool is_euler_cromer{};
    
    Position position{};

    Velocity velocity{polar(constants::start_speed, constants::start_angle)};

    void step(units::physical::Time auto const time_step) 
    {
        auto const update_position = [=, this] { position += velocity*time_step; };
        
        if (!is_euler_cromer) {
            update_position();
        }

        if (is_damped) {
            auto const damping_force = constants::air_resistance*velocity.length_squared();
            auto const damping_acceleration = -velocity.normalized()*(damping_force/constants::ball_mass);
            velocity += (damping_acceleration + constants::g)*time_step;
        }
        else {
            velocity += constants::g*time_step;
        }

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
    return matplot::plot(x, y)->line_width(1.f);
}

void plot_all() 
{
    matplot::colororder(std::vector<std::vector<float>>{
        {1.f, 0.4f, 0.4f}, {0.4f, 0.f, 0.f},
        {0.4f, 1.f, 0.4f}, {0.f, 0.4f, 0.f}
    });

    for (auto const time_step : constants::time_steps) 
    {
        plot_simulation({.is_damped = false, .is_euler_cromer = false}, time_step);
        plot_simulation({.is_damped = false, .is_euler_cromer = true}, time_step);
        plot_simulation({.is_damped = true, .is_euler_cromer = false}, time_step);
        plot_simulation({.is_damped = true, .is_euler_cromer = true}, time_step);
    }

    auto const legend = matplot::legend({"Odämpad", "Odämpad, Euler-Cromer", "Dämpad", "Dämpad, Euler-Cromer"});
    legend->box(false);
    legend->location(matplot::legend::general_alignment::bottomright);
}

//------------------------------------------

int main() 
{
    matplot::gcf()->reactive_mode(false);
    matplot::gcf()->size(1600, 700);
    matplot::hold(true);
    matplot::grid(true);
    
    matplot::xlabel("X-position [meter]");
    matplot::ylabel("Y-position [meter]");

    plot_all();

    matplot::axis(matplot::tight);

    matplot::save("results/pingisboll.svg");
    matplot::save("results/pingisboll.jpeg");
}
