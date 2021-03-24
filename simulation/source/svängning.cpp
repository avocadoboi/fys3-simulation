// Björn Sundin, 2021

#include <common.hpp>

//------------------------------------------

namespace constants {

constexpr auto mass = 0.3q_kg;
constexpr auto spring_constant = 15.q_N/1.q_m;
constexpr auto damping_ratio = 2.l;

auto const critical_damping = 2.*units::sqrt(mass*spring_constant);
auto const damping_coefficient = damping_ratio*critical_damping;

constexpr auto start_pos = 1.q_m;
constexpr auto start_velocity = 0.q_m_per_s;

constexpr auto time_steps = std::array{0.1q_ms, 1.q_ms, 5.q_ms, 10.q_ms};
constexpr auto time = 1.5q_s;

} // namespace constants

//------------------------------------------

struct Simulation 
{
    bool is_euler_cromer{};
    PositionScalar y{constants::start_pos};
    VelocityScalar v{constants::start_velocity};
    AccelerationScalar a{calculate_acceleration()};
    si::time<si::second, long double> t;

    AccelerationScalar calculate_acceleration() {
        return -(y*constants::spring_constant + v*constants::damping_coefficient)/constants::mass;
    }

    void step(units::physical::Time auto const time_step) 
    {
        auto const update_position = [=, this] { y += units::quantity_cast<si::metre>(v*time_step); };

        if (!is_euler_cromer) {
            update_position();
        }

        v += units::quantity_cast<si::metre_per_second>(a*time_step);
        
        if (is_euler_cromer) {
            update_position();
        }
        
        a = calculate_acceleration();

        t += units::quantity_cast<si::second>(time_step);
    }
};

//------------------------------------------

void plot_simulation(Simulation simulation, units::physical::Time auto const time_step) 
{
    auto const iterations = calculate_iterations(time_step, constants::time);
    auto const new_data = [=]{ return std::vector<long double>(iterations); };

    auto time_data = new_data();
    auto position_data = new_data();
    auto velocity_data = new_data();
    auto acceleration_data = new_data();

    for (auto const i : range(iterations)) 
    {
        time_data[i] = simulation.t.count();
        position_data[i] = simulation.y.count();
        velocity_data[i] = simulation.v.count();
        acceleration_data[i] = simulation.a.count();

        simulation.step(time_step);
    }

    matplot::plot(time_data, position_data)->line_width(1.f);
    matplot::plot(time_data, velocity_data)->line_width(1.f);
    matplot::plot(time_data, acceleration_data)->line_width(1.f).use_y2(true);
    matplot::gca()->y2_axis().color(matplot::color::blue).touch();
}

void plot_all() 
{
    matplot::colororder(std::vector<std::vector<float>>{
        {1.f, 0.4f, 0.4f}, {0.4f, 1.f, 0.4f}, {0.4f, 0.4f, 1.f},
        {0.4f, 0.f, 0.f}, {0.f, 0.4f, 0.f}, {0.f, 0.f, 0.4f}
    });
    matplot::gca()->y_axis().color({1.f, 0.f, 0.f});
    matplot::gca()->y2_axis().color({0.f, 0.f, 1.f});
    
    for (auto const time_step : constants::time_steps) 
    {
        plot_simulation({.is_euler_cromer = false}, time_step);
        plot_simulation({.is_euler_cromer = true}, time_step);
    }

    auto const legend = matplot::legend({
        "Position [m]", "Fart [m/s]", "Acceleration [m/s^2]",
        "Position, Euler-Cromer [m]", "Fart, Euler-Cromer [m/s]", "Acceleration, Euler-Cromer [m/s^2]",
    });
    legend->box(false);
    legend->location(matplot::legend::general_alignment::bottomright);
}

std::string get_file_name() {
    auto name = fmt::format("results/svängning_dämpningsratio_{}", constants::damping_ratio);
    std::ranges::replace(name, '.', '_');
    return name;
}

//------------------------------------------

int main() 
{
    matplot::gcf()->reactive_mode(false);
    matplot::gcf()->size(1600, 700);
    matplot::hold(true);
    // matplot::grid(true);

    matplot::xlabel("Tid [s]");
    matplot::ylabel("Position [m] och fart [m/s]");
    matplot::y2label("Acceleration [m/s^2]");
    
    plot_all();

    auto const name = get_file_name();
    matplot::save(name, "jpeg");
    matplot::save(name, "svg");
}
