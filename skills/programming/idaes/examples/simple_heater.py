"""
Simple heater example using IDAES.

This example demonstrates:
- Creating a basic flowsheet
- Adding a property package (IAPWS95 for water/steam)
- Adding a heater unit model
- Specifying inputs
- Initializing the model
- Solving
- Displaying results

Problem: Heat water from liquid state to saturated steam
- Inlet: 100 mol/s liquid water at 25°C, 1 atm
- Heat duty: 500 kW
- Find outlet conditions
"""

from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom


def main():
    # Create Pyomo model and IDAES flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Add property package for water/steam (IAPWS95)
    m.fs.properties = iapws95.Iapws95ParameterBlock()

    # Add heater unit model
    m.fs.heater = Heater(property_package=m.fs.properties)

    # Specify inlet conditions
    # Flow: 100 mol/s
    m.fs.heater.inlet.flow_mol[0].fix(100)  # mol/s

    # Pressure: 1 atm = 101325 Pa
    m.fs.heater.inlet.pressure[0].fix(101325)  # Pa

    # Temperature: 25°C = 298.15 K
    # IAPWS uses (flow, enthalpy, pressure) state vars
    # Need to calculate enthalpy for liquid water at 25°C, 1 atm
    # Approximate: h ≈ 104.8 kJ/mol = 104800 J/mol
    m.fs.heater.inlet.enth_mol[0].fix(104800)  # J/mol

    # Specify heat duty: 500 kW = 500000 W
    m.fs.heater.heat_duty[0].fix(500000)  # W

    # Check degrees of freedom (should be 0)
    dof = degrees_of_freedom(m)
    print(f"Degrees of Freedom: {dof}")
    assert dof == 0, "Model not fully specified"

    # Initialize the heater
    print("\nInitializing heater...")
    m.fs.heater.initialize()
    print("Initialization complete.")

    # Solve the model
    print("\nSolving model...")
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # Check solver status
    from pyomo.opt import TerminationCondition

    if results.solver.termination_condition == TerminationCondition.optimal:
        print("\n" + "=" * 60)
        print("Solve successful!")
        print("=" * 60)
    else:
        print(f"\nWarning: Solver returned status: {results.solver.termination_condition}")

    # Display results
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)

    # Inlet conditions
    print("\nInlet Conditions:")
    print(f"  Flow: {value(m.fs.heater.inlet.flow_mol[0]):.2f} mol/s")
    inlet_p = value(m.fs.heater.inlet.pressure[0])
    print(f"  Pressure: {inlet_p:.0f} Pa ({inlet_p / 101325:.2f} atm)")
    print(f"  Enthalpy: {value(m.fs.heater.inlet.enth_mol[0]):.0f} J/mol")

    # Outlet conditions
    print("\nOutlet Conditions:")
    print(f"  Flow: {value(m.fs.heater.outlet.flow_mol[0]):.2f} mol/s")
    outlet_p = value(m.fs.heater.outlet.pressure[0])
    print(f"  Pressure: {outlet_p:.0f} Pa ({outlet_p / 101325:.2f} atm)")
    print(f"  Enthalpy: {value(m.fs.heater.outlet.enth_mol[0]):.0f} J/mol")

    # Heat duty
    print("\nHeat Duty:")
    print(f"  {value(m.fs.heater.heat_duty[0]) / 1000:.1f} kW")

    # Energy balance check
    h_in = value(m.fs.heater.inlet.enth_mol[0])
    h_out = value(m.fs.heater.outlet.enth_mol[0])
    flow = value(m.fs.heater.inlet.flow_mol[0])
    energy_change = (h_out - h_in) * flow
    print(f"\nEnergy Balance Check:")
    print(f"  ΔH = {(h_out - h_in) / 1000:.1f} kJ/mol")
    print(f"  Energy added = {energy_change / 1000:.1f} kW")
    print(f"  Heat duty = {value(m.fs.heater.heat_duty[0]) / 1000:.1f} kW")
    print(f"  Difference = {abs(energy_change - value(m.fs.heater.heat_duty[0])):.2e} W")

    # Full unit report
    print("\n" + "=" * 60)
    print("DETAILED UNIT REPORT")
    print("=" * 60)
    m.fs.heater.report()

    return m


if __name__ == "__main__":
    model = main()
