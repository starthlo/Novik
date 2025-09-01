#!/usr/bin/env python3
"""
Test script for user profile endpoints
"""
import requests
import json

BASE_URL = "http://localhost:8000/api"

def test_user_endpoints():
    print("Testing User Profile Endpoints\n")
    print("=" * 50)
    
    # Step 1: Login to get access token
    print("\n1. Testing Login...")
    login_data = {
        "email": "test@example.com",  # You'll need to use a valid user
        "password": "testpassword"
    }
    
    try:
        response = requests.post(f"{BASE_URL}/auth/login/", json=login_data)
        
        if response.status_code == 200:
            data = response.json()
            access_token = data.get("accessToken")
            user_info = data.get("user")
            print(f"✓ Login successful!")
            print(f"  User: {user_info.get('email')}")
            print(f"  Token: {access_token[:20]}...")
            
            # Step 2: Test Get Profile endpoint
            print("\n2. Testing Get Profile...")
            headers = {"Authorization": f"Bearer {access_token}"}
            profile_response = requests.get(f"{BASE_URL}/user/profile/", headers=headers)
            
            if profile_response.status_code == 200:
                profile = profile_response.json()
                print(f"✓ Profile retrieved successfully!")
                print(f"  Username: {profile.get('username')}")
                print(f"  Email: {profile.get('email')}")
                print(f"  Name: {profile.get('first_name')} {profile.get('last_name')}")
            else:
                print(f"✗ Failed to get profile: {profile_response.status_code}")
                print(f"  Response: {profile_response.text}")
            
            # Step 3: Test Update Profile endpoint
            print("\n3. Testing Update Profile...")
            update_data = {
                "first_name": "Test",
                "last_name": "User",
                "occupation": "Dentist",
                "city": "New York"
            }
            
            update_response = requests.patch(
                f"{BASE_URL}/user/profile/update/", 
                json=update_data, 
                headers=headers
            )
            
            if update_response.status_code == 200:
                updated_profile = update_response.json()
                print(f"✓ Profile updated successfully!")
                print(f"  Updated Name: {updated_profile.get('first_name')} {updated_profile.get('last_name')}")
                print(f"  Occupation: {updated_profile.get('occupation')}")
            else:
                print(f"✗ Failed to update profile: {update_response.status_code}")
                print(f"  Response: {update_response.text}")
                
        else:
            print(f"✗ Login failed: {response.status_code}")
            print(f"  Response: {response.text}")
            print("\nNote: You need to create a test user first or use valid credentials")
            
    except requests.exceptions.ConnectionError:
        print("✗ Could not connect to server. Make sure the Django server is running on port 8000")
    except Exception as e:
        print(f"✗ Error: {e}")

if __name__ == "__main__":
    test_user_endpoints()
    print("\n" + "=" * 50)
    print("Test completed!")
