"""
Flash separation example using IDAES.

This example demonstrates:
- Creating a flowsheet with multiple units
- Connecting units with Arcs
- Using flash separator for vapor-liquid separation
- Sequential initialization
- Analyzing separation performance

Problem: Separate a heated water stream using flash separation
- Feed: 100 mol/s water at high enthalpy
- Flash at lower pressure to create vapor and liquid products
- Calculate vapor and liquid flow rates and compositions
"""

from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater, Flash
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import propagate_state
from pyomo.network import Arc, expand_arcs


def main():
    # Create model and flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Property package
    m.fs.properties = iapws95.Iapws95ParameterBlock()

    # Add unit models
    # Heater to bring water to desired conditions
    m.fs.heater = Heater(property_package=m.fs.properties)

    # Flash separator for vapor-liquid separation
    m.fs.flash = Flash(
        property_package=m.fs.properties, has_heat_transfer=True, has_pressure_change=True
    )

    # Connect heater outlet to flash inlet
    m.fs.stream_01 = Arc(source=m.fs.heater.outlet, destination=m.fs.flash.inlet)

    # Expand arcs to create connection constraints
    expand_arcs(m)

    # Specify heater inlet (subcooled liquid water)
    m.fs.heater.inlet.flow_mol[0].fix(100)  # mol/s
    m.fs.heater.inlet.pressure[0].fix(500000)  # 5 bar = 500000 Pa
    m.fs.heater.inlet.enth_mol[0].fix(104800)  # ~25°C liquid water, J/mol

    # Heat to higher enthalpy
    m.fs.heater.heat_duty[0].fix(1000000)  # 1 MW

    # Flash specifications
    # Lower pressure to promote flashing
    m.fs.flash.deltaP[0].fix(-300000)  # Drop to ~2 bar

    # Adiabatic flash (no heat transfer)
    m.fs.flash.heat_duty[0].fix(0)

    # Check DOF
    dof = degrees_of_freedom(m)
    print(f"Degrees of Freedom: {dof}")
    assert dof == 0, "Model not fully specified"

    # Initialize sequentially
    print("\nInitializing heater...")
    m.fs.heater.initialize()

    print("Propagating state to flash...")
    propagate_state(arc=m.fs.stream_01)

    print("Initializing flash...")
    m.fs.flash.initialize()

    print("Initialization complete.\n")

    # Solve the model
    print("Solving flowsheet...")
    solver = get_solver()
    results = solver.solve(m, tee=False)

    # Check solver status
    from pyomo.opt import TerminationCondition

    if results.solver.termination_condition == TerminationCondition.optimal:
        print("Solve successful!\n")
    else:
        print(f"Warning: Solver status: {results.solver.termination_condition}\n")

    # Display results
    print("=" * 60)
    print("FLASH SEPARATION RESULTS")
    print("=" * 60)

    # Feed conditions (flash inlet)
    print("\nFeed to Flash:")
    print(f"  Flow: {value(m.fs.flash.inlet.flow_mol[0]):.2f} mol/s")
    print(f"  Pressure: {value(m.fs.flash.inlet.pressure[0]) / 1000:.1f} kPa")
    print(f"  Enthalpy: {value(m.fs.flash.inlet.enth_mol[0]) / 1000:.1f} kJ/mol")

    # Vapor outlet
    print("\nVapor Product:")
    flow_vap = value(m.fs.flash.vap_outlet.flow_mol[0])
    print(f"  Flow: {flow_vap:.2f} mol/s")
    print(f"  Pressure: {value(m.fs.flash.vap_outlet.pressure[0]) / 1000:.1f} kPa")
    print(f"  Enthalpy: {value(m.fs.flash.vap_outlet.enth_mol[0]) / 1000:.1f} kJ/mol")

    # Liquid outlet
    print("\nLiquid Product:")
    flow_liq = value(m.fs.flash.liq_outlet.flow_mol[0])
    print(f"  Flow: {flow_liq:.2f} mol/s")
    print(f"  Pressure: {value(m.fs.flash.liq_outlet.pressure[0]) / 1000:.1f} kPa")
    print(f"  Enthalpy: {value(m.fs.flash.liq_outlet.enth_mol[0]) / 1000:.1f} kJ/mol")

    # Separation performance
    print("\nSeparation Performance:")
    feed_flow = value(m.fs.flash.inlet.flow_mol[0])
    vapor_fraction = (flow_vap / feed_flow) * 100
    liquid_fraction = (flow_liq / feed_flow) * 100
    print(f"  Vapor fraction: {vapor_fraction:.1f}%")
    print(f"  Liquid fraction: {liquid_fraction:.1f}%")

    # Mass balance check
    total_out = flow_vap + flow_liq
    print(f"\nMass Balance:")
    print(f"  Feed flow: {feed_flow:.2f} mol/s")
    print(f"  Total outlet: {total_out:.2f} mol/s")
    print(f"  Difference: {abs(feed_flow - total_out):.2e} mol/s")

    # Energy balance check
    h_in = value(m.fs.flash.inlet.enth_mol[0])
    h_vap = value(m.fs.flash.vap_outlet.enth_mol[0])
    h_liq = value(m.fs.flash.liq_outlet.enth_mol[0])
    q_flash = value(m.fs.flash.heat_duty[0])

    energy_in = feed_flow * h_in + q_flash
    energy_out = flow_vap * h_vap + flow_liq * h_liq

    print(f"\nEnergy Balance:")
    print(f"  Energy in: {energy_in / 1000:.1f} kW")
    print(f"  Energy out: {energy_out / 1000:.1f} kW")
    print(f"  Difference: {abs(energy_in - energy_out):.2e} W")

    # Full flash report
    print("\n" + "=" * 60)
    print("DETAILED FLASH REPORT")
    print("=" * 60)
    m.fs.flash.report()

    return m


if __name__ == "__main__":
    model = main()
