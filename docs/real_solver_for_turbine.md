For a **basic CFD solver coupled to a turbine**, I’d avoid trying to make CFD do everything. The most useful architecture is a **reduced-order turbine model coupled to CFD**.

For a wind or water turbine, the core coupling would be:

**CFD → forces/torque → rotor dynamics → rotor state → CFD**

### Minimal useful stack

| Component | Purpose | Complexity |
|---|---|---:|
| **CFD solver** | Velocity, pressure, turbulence/viscosity field | ★★ |
| **Blade geometry** | Airfoil/hydrofoil sections + 3D blade | ★★ |
| **Force model** | Pressure + viscous forces on blades | ★★ |
| **Rotor model** | Angular velocity, acceleration, torque | ★ |
| **Generator/load model** | Converts torque → electrical/mechanical load | ★ |
| **Control model** | Pitch, torque, RPM control | ★★ |
| **BEM / blade-element model** | Fast approximation for comparison/initialization | ★ |
| **Motion/mesh model** | Rotor rotation and possibly moving/deforming mesh | ★★★ |

The **minimum viable turbine simulation** could actually be:

```text
                 ┌──────────────┐
                 │ Blade/Geometry│
                 └──────┬───────┘
                        │
                        ▼
┌──────────┐     ┌──────────────┐
│ CFD Field│────►│ Blade Forces │
└────▲─────┘     └──────┬───────┘
     │                  │
     │                  ▼
     │           ┌──────────────┐
     │           │ Torque Model │
     │           └──────┬───────┘
     │                  │
     │                  ▼
     │           ┌──────────────┐
     └───────────│ Rotor Dynamics│
                 └──────┬───────┘
                        │
                        ▼
                  RPM / position
                        │
                        └──────► CFD
```


