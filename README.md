# APSIS

**APSIS** is an interpreted domain-specific programming language designed specifically for scientific computing, celestial mechanics, orbital dynamics, and numerical analysis.

___

## Features

- **Built-in Astrodynamics Engine**: First-class support for state vector propagation, orbital element conversion (Keplerian, Cartesian), and trajectory calculations.
- **Scientific Unit Aware**: Built-in handling for units and dimensional analysis.
- **Linear Algebra & Numerics**: Fast vector/matrix operations and numerical integration methods.
- **Python-like CLI**: Execute scripts seamlessly from anywhere on your machine using `apsis <script.jdx>`.

___

## Installation

### Windows (Quick Installer)

1. Clone or download this repository to your computer:
   ```bash
   git clone https://github.com/Jits-Doomen/Apsis.git
   cd Apsis
   ```

2. Run the automated installer:
   - Double-click **`installer.bat`**, or run it from your command prompt:
     ```cmd
     .\installer.bat
     ```

3. **Restart your terminal** or reload your environment variables:
   ```powershell
   # In PowerShell (or simply close and reopen your terminal)
   $env:Path = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
   ```

4. Verify installation:
   ```bash
   apsis
   ```

___

## Building from Source (Manual)

### Prerequisites
- **CMake** (v3.16 or higher)
- **C++20 Compliant Compiler** (MSVC, GCC, or Clang)

### Build Steps

```bash
# 1. Clone repository
git clone https://github.com/Jits-Doomen/Apsis.git
cd Apsis

# 2. Configure build system
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release

# 3. Compile project
cmake --build build --config Release
```

The compiled binary will be generated at `build/Release/apsis.exe`.

___

## Usage

Run any `.jdx` script directly from any directory:

```bash
apsis path/to/script.jdx
```

### Example Script (`orbit_demo.jdx`)

```apsis
// Calculate circular orbital velocity at 400 km (~248,548 miles) LEO
var R_earth = 6378.137;
var altitude = 400.0;
var mu_earth = 398600.4418;

var r = R_earth + altitude;
var v = sqrt(mu_earth / r);

print("Orbital Radius: " + r + " km");
print("Circular Velocity: " + v + " km/s");
```

Run it:
```bash
apsis orbit_demo.jdx
```

___

## Repository Structure

```text
Apsis/
├── include/              # Public C++ header files (.hpp)
│   └── apsis/            # Language interpreter & core math headers
├── src/                  # Source implementation files (.cpp)
├── packages/             # Physics, astrodynamics, and ephemerides modules
├── examples/             # Sample .jdx scripts
├── CMakeLists.txt        # Primary CMake build configuration
├── installer.bat         # One-click Windows installer script
└── README.md             # Project documentation
```

___

## Contributing

Contributions are welcome! If you find a bug or want to add a feature:
1. Fork the repository.
2. Create your feature branch (`git checkout -b feature/yourchange`).
3. Commit your changes (`git commit -m 'Add new (insert your change here)'`).
4. Push to the branch (`git push origin feature/yourchange`).
5. Open a Pull Request.

___

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
