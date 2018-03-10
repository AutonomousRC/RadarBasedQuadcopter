----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date: 2018/03/09 14:48:58
-- Design Name: 
-- Module Name: servo_pwm - Behavioral
-- Project Name: 
-- Target Devices: 
-- Tool Versions: 
-- Description: 
-- 
-- Dependencies: 
-- 
-- Revision:
-- Revision 0.01 - File Created
-- Additional Comments:
-- 
----------------------------------------------------------------------------------


library IEEE;                                       --라이브러리, 패키지선언
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity servo_pwm is                                 -- entity 선언
    PORT (
        clk   : IN  STD_LOGIC;
        reset : IN  STD_LOGIC;
        pos   : IN  STD_LOGIC_VECTOR(7 downto 0);
        servo : OUT STD_LOGIC
    );
end servo_pwm;

architecture Behavioral of servo_pwm is             
                                                     -- Counter, from 0 to 3599. 20ms를 만들기 위해서 180x20 = 3600이 필요
    signal cnt : unsigned(11 downto 0);
                                                     -- Temporal signal used to generate the PWM pulse.
    signal pwmi: unsigned(9 downto 0);
begin
                                                     -- 최소 duty cycle이 1ms가 되어야 하므로 180을 더해줌.
    pwmi <= "00" & unsigned(pos) + 180;
                                                     -- Counter process, from 0 to 3599.
    counter: process (reset, clk) begin
        if (reset = '1') then                        -- reset 설정
            cnt <= (others => '0');
        elsif rising_edge(clk) then                  -- cnt=3599에 도달하면 cnt를 0으로!
            if (cnt = 3599) then
                cnt <= (others => '0');
            else
                cnt <= cnt + 1;                      -- cnt < 3599 이면 +1
            end if;
        end if;
    end process;
                                                     -- pwm 출력
    servo <= '1' when (cnt < pwmi) else '0';
end Behavioral;
