#pragma once

#include <cstdint>

static constexpr uint64_t FILL_NUMBER_MASK = 0x000000000000FFFF;
static constexpr uint64_t ACCELERATOR_MODE_MASK = 0x00000000001F0000;
static constexpr uint64_t BEAM_MODE_MASK = 0x0000000003E00000;
static constexpr uint64_t FAST_FILTER_MASK = 0x00000001FC000000;
static constexpr uint64_t ADVANCED_FILTER_MASK = 0x00001FFE00000000;

static constexpr uint64_t FAST_FILTER_SCIFI = 0x1;
static constexpr uint64_t FAST_FILTER_SCIFI_TOTAL = 0x2;
static constexpr uint64_t FAST_FILTER_US = 0x4;
static constexpr uint64_t FAST_FILTER_US_TOTAL = 0x8;
static constexpr uint64_t FAST_FILTER_DS = 0x10;
static constexpr uint64_t FAST_FILTER_DS_TOTAL = 0x20;
static constexpr uint64_t FAST_FILTER_VETO_TOTAL = 0x40;

static constexpr uint64_t ADVANCED_FILTER_SCIFI_PLANES = 0x1;
static constexpr uint64_t ADVANCED_FILTER_SCIFI_HITS = 0x2;
static constexpr uint64_t ADVANCED_FILTER_US_PLANES = 0x4;
static constexpr uint64_t ADVANCED_FILTER_US_HITS = 0x8;
static constexpr uint64_t ADVANCED_FILTER_DSH_PLANES = 0x10;
static constexpr uint64_t ADVANCED_FILTER_DSH_HITS = 0x20;
static constexpr uint64_t ADVANCED_FILTER_DSV_PLANES = 0x40;
static constexpr uint64_t ADVANCED_FILTER_DSV_HITS = 0x80;
static constexpr uint64_t ADVANCED_FILTER_DS_PLANES = 0x100;
static constexpr uint64_t ADVANCED_FILTER_VETO_PLANES = 0x200;
static constexpr uint64_t ADVANCED_FILTER_VETO_HITS = 0x400;
static constexpr uint64_t ADVANCED_FILTER_GLOBAL_PLANES = 0x800;

enum class LhcAcceleratorMode {
  Null,
  Shutdown,
  Cooldown,
  MachineCheckout,
  Access,
  MachineTest,
  Calibration,
  WarmUp,
  Recovery,
  SectorDependent,
  BeamSetup,
  ProtonPhysics,
  IonPhysics,
  TotemPhysics,
  SpecialOpticsPhysics,
  ProtonNucleusPhysics,
  MachineDevelopment,
  Unknown = 31
};

enum class LhcBeamMode {
  Null,
  Setup,
  Abort,
  InjectionProbeBeam,
  InjectionSetupBeam,
  InjectionPhysicsBeam,
  PrepareRamp,
  Ramp,
  FlatTop,
  Squeeze,
  Adjust,
  StableBeams,
  UnstableBeams,
  BeamDump,
  RampDown,
  Cycling,
  Recovery,
  InjectAndDump,
  CirculateAndDump,
  NoBeam,
  Unknown = 31
};
