----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date: 2018/03/09 15:30:08
-- Design Name: 
-- Module Name: servo_pwm_180kHZ - Behavioral
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


library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity servo_pwm_clk180kHz is                           -- entity
    PORT(
        clk  : IN  STD_LOGIC;
        reset: IN  STD_LOGIC;
        sw : in std_logic_vector(3 downto 0);
        
        servo: OUT STD_LOGIC
    );
end servo_pwm_clk180kHz;

architecture Behavioral of servo_pwm_clk180kHz is       --architecture
    component swt 
--  Port ( );
    port(   
        sw : in std_logic_vector(3 downto 0);
        pos : out std_logic_vector(7 downto 0)
        );
        
    end component;
    
    COMPONENT clk_180kHz                                --clk_180kHz 의 인스턴스 생성
        PORT(
            clk    : in  STD_LOGIC;
            reset  : in  STD_LOGIC;
            clk_out: out STD_LOGIC
        );
    END COMPONENT;
    
    COMPONENT servo_pwm                                 --servo_pwm의 인스턴스 생성
        PORT (
            clk   : IN  STD_LOGIC;
            reset : IN  STD_LOGIC;
            pos   : IN  STD_LOGIC_VECTOR(7 downto 0);
            servo : OUT STD_LOGIC
        );
    END COMPONENT;
    
    signal clk_out : STD_LOGIC := '0';
    signal pos  : STD_LOGIC_VECTOR(7 downto 0);
begin
    swt_in : swt port map(
        sw=>sw, pos => pos
        );
        

    clk180kHz_map: clk_180kHz PORT MAP(                 --인스턴스 이름 : 하위레벨 entity 이름
        clk=>clk, reset=>reset, clk_out=>clk_out        --포트 연결
    );
    
    servo_pwm_map: servo_pwm PORT MAP(                      --인스턴스 이름 : 하위레벨 entity 이름
        clk=>clk_out, reset=>reset, pos=>pos, servo=>servo  --포트 연결
    );
end Behavioral;
