# 5G Multi-layer NOMA Transmission Simulation System for Mobile user

## Introduction

This project provides a simulation system for multi-layer Non-Orthogonal Multiple Access (NOMA) transmission within the 5G downlink physical layer, in line with the 3GPP TR 38.901 specification. Our system is designed to process user data concurrently, leveraging the DL-SCH and PDSCH encoding schemes, and multiplexing the data using a NOMA-specific method dictated by predefined power allocation (PA) parameters. 

The system integrates complex functionalities such as MIMO precoding, OFDM modulation, and the use of a single channel for communicating with multiple receivers. It also features a critical buffer at the receiver end for storing composite NOMA signals, which is vital for successive interference cancellation (SIC).

## System Model

![Simulation Model](https://github.com/iamsip/NOMA5G_MOBILITY/blob/main/5GwithPA.jpg)

The system model, as shown in the figure above, demonstrates the flow of data from transmission to reception. At the core of the transmission process is the NOMA multiplexer, which handles data according to PA parameters. The transmitted NOMA signal undergoes MIMO precoding and OFDM modulation before it is sent to the receivers.

## Mobility Analysis

The simulation system also factors in mobility analysis for varying transmission environments, calculating 2D and 3D distances between the base station (BS) and the user. This analysis further extends to path loss calculation, reliant on the derived distances, channel frequency, and the specific transmission environment, be it urban or rural.

## Power Allocation and Performance Optimization
The system employs bit error rate (BER) feedback from receivers, allowing the transmitter to adjust PA in an iterative manner across layers for future transmissions. This approach aims to enhance system throughput, meet BER prerequisites, and maintain dependable NOMA connections, even under unstable channel conditions.

## Goals

- Maximize system throughput.
- Ensure adherence to BER standards.
- Maintain robust NOMA connections in varying channel conditions.

## Getting Started

Matlab 5G Toolbox is required

## References

1. [3GPP TR 38.901 - "Study on channel model for frequencies from 0.5 to 100 GHz"](https://www.3gpp.org/DynaReport/38901.htm)
