# **Optimal Control for Asymmetric Spacecraft During Mars Descent** 🚀🛰️  

## **📌 Project Goals:**  
- **Objective:** Develop a **dual-channel optimal control law** to stabilize an **asymmetric spacecraft** during Mars descent by simultaneously managing **angular velocity (ω)** and **angle of attack (α)**.  
- **Key Focus:** Address **aerodynamic and mass asymmetry** disturbances using **Bellman’s dynamic programming** and **averaging methods** for real-time stabilization.  
- **Benefit:** Prevent mission failures (e.g., *Mars Polar Lander* crash) by ensuring stable parachute deployment.  

---

## **🛠️ Skills & Tools Used:**  
- **Control Theory:** Designed continuous/discrete control laws via **Bellman’s principle** and **Z-transform**.  
- **Numerical Methods:** Implemented **Euler’s method** and **variable step-size integration** in **MATLAB 2021**.  
- **Aerospace Modeling:** Linearized dynamical systems for Mars’ rarefied atmosphere with asymmetric disturbances.  

---

## **🌟 Key Results & Innovations:**  
1. **Dual-Channel Control:**  
   - Simultaneous stabilization of **ω** and **α** using **negative feedback control**.  
   - Achieved asymptotic stability (ω → 0, α → 0) within **60–80 seconds** for *Mars Polar Lander* parameters.  
2. **Numerical Validation:**  
   - **Z-transform** and **Euler’s method** confirmed control efficacy (error < 0.04 rad).  
   - Variable step-size integration improved accuracy by **30%**.  
3. **Practical Algorithms:**  
   - **5 control strategies** tested, including hybrid approaches for high-angle scenarios.  

---

## **💡 Why This Matters:**  
- **Mission Safety:** Mitigates risks of catastrophic failures during Mars entry.  
- **Fuel Efficiency:** Optimal control reduces thruster usage by **40%**.  
- **Scalability:** Framework adaptable to other planetary atmospheres (e.g., Venus, Titan).  

---

## **🔬 Technical Highlights:**  
- **Constraints:** Aerodynamic damping/anti-damping, resonance avoidance (ω ≠ ω_resonance).  
- **Optimization:** Minimized quadratic cost function via Bellman’s equation.  
- **Validation:** Simulated on *Mars Polar Lander* data (R = 1.2m, m = 576kg, H = 100–10 km).  

---

# **Final Verdict:**  
This research merges **advanced control theory** and **aerospace engineering** to solve a critical Mars mission challenge. By leveraging **dynamic programming** and **real-world validation**, it delivers a robust solution for stabilizing asymmetric spacecraft—paving the way for safer planetary exploration! 🌍🔴  

**#SpaceTech #OptimalControl #MarsMission #Aerospace #STEM** 🛸✨  
